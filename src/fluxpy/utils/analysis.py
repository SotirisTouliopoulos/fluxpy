
import numpy as np
import pandas as pd
import math
import networkx as nx
from .utils import NestedDataFrameType

def get_nutrients_gradient(model, nutrients=None, upper_bound=None, step=None) -> NestedDataFrameType:
    """
    Assuming a single medium comound changes at a time, it returns FVA findings over the chaning envs.
    The model.medium will be used, i.e. first set the medium you want to investigate to your model.
    Returns a df with as many rows as the model.medium and columns as many as their quotient by the step.
    Each cell of this df includes the fva result for the medium with the corresponding alteration.
    For example, df.loc[0,0] corresponds to the medium where the first compound has the lowest value possible.

    model -- a cobra model
    nutrients -- list of reaction ids for which a gradient will be calculated
    upper_bound -- maximum a flux can get
    step -- increase of upper_bound in each iteration of the gradient
    """
    from cobra.flux_analysis import flux_variability_analysis
    if upper_bound is None:
        upper_bound = max(model.medium.values())
    if step is None:
        num_of_chunks = 5
        step = math.floor(upper_bound / num_of_chunks)
    else:
        num_of_chunks = math.floor(upper_bound / step)

    modulo = (upper_bound % num_of_chunks) + 1

    if nutrients is None:
        nutrients = model.medium.copy()
        num_of_medium_compounds = len(nutrients)
    else:
        num_of_medium_compounds = len(nutrients)

    # [NOTE] If a column is added at a later point using the pd.at() it's gonna be nan
    df = pd.DataFrame(index=range(num_of_medium_compounds), columns=range(num_of_chunks + 1))

    initial_model = model.copy()
    x_coord = 0

    for compound in nutrients:
        print(compound, str(x_coord) + "/" + str(num_of_medium_compounds))
        model_in = initial_model.copy()
        j_coord = 0
        for j in np.arange(0, upper_bound+modulo, step):
            if j > upper_bound:
                j = upper_bound
            limited_medium = model_in.medium
            limited_medium[compound] = j
            model_in.medium = limited_medium
            try:
                fva = flux_variability_analysis(model_in)
            except:
                print("Infeasible for ", compound, "with upper bound of ", str(j))
                fva = pd.DataFrame(data=[[0]*len(model.reactions)] * 2).T
                fva.columns = ["minimum", "maximum"]
                fva.index =  [x.id for x in model.reactions]
            df.at[x_coord, j_coord] = fva
            j_coord += 1
        x_coord += 1

    df.index = nutrients
    columns = []
    for i in range(num_of_chunks + 1):
        columns.append("ub_" + str(i*step))
    df.columns = columns

    return df


# %% Util functions: flux coupling analysis
def parse_qfca(qfca_output, model=None, remove_exchange_routes=True, exclude_biomass=True, format="csv"):
    """
    Parses QFCA (Quantitative Fatty Acid Composition Analysis) output data into a networkx graph representation.

    Args:
        qfca_output (str): Path to the QFCA output file in CSV or Excel format.
        model (object, optional): Metabolic model object. Default is None.
        remove_exchange_routes (bool, optional): Whether to remove exchange reactions from the graph. Default is True.
        exclude_biomass (bool, optional): Whether to exclude biomass reaction from the graph. Default is True.
        format (str, optional): Format of the QFCA output file. Can be 'csv' or 'xlsx'. Default is 'csv'.

    Returns:
        nx.Graph: A networkx Graph representing the interactions between metabolic reactions based on QFCA data.
    """
    if format=="csv":
        df = pd.read_csv(qfca_output, index_col=0)
    elif format=="xlsx":
        df = pd.read_excel(qfca_output, index_col=0)

    # Create a graphimport networkx as nx
    G = nx.Graph()

    # Add nodes
    G.add_nodes_from(df.index)

    # Check if model is provided and get biomass_reaction if possible
    biomass_reaction = None
    if model is not None and model.objective is not None:
        objective = str(model.summary()._objective)
        biomass_reaction = objective.split(" ")[1]

    """
    1 - fully coupled reactions - Black // Grey // dark green (#004D40)
    2 - partially coupled reactions - Blue // SteelBlue // Cyan  (#1E88E5)
    3 - reaction i is directionally coupled to reaction j - Red // DarkSalmon // magenta (#D81B60)
    4 - reaction j is directionally coupled to reaction i - Green // LightSeaGreen // dark yellow  (#FFC107)

    [NOTE] we now use dark yellow for both 3 and 4 and the magenta as EDGE_SECTED_PAINT in the style json file.
    """
    pallette = {1: 'black',
                2: 'blue',
                3: 'gree',
                4: 'green'
    }
    parsed_pairs = set()
    counter = 0
    for source in df.index:
        if model is not None and exclude_biomass and biomass_reaction is not None:
            if source == biomass_reaction:
                continue
        for target in df.index:
            if model is not None and exclude_biomass and biomass_reaction is not None:
                if target == biomass_reaction:
                    continue
            # Check if pair is already parsed
            if (source, target) in parsed_pairs or (target, source) in parsed_pairs:
                continue
            # Add pair in parsed and
            parsed_pairs.add((source, target))
            value = df.loc[source, target]
            # Ignore cases where i = j or value is 0
            if source != target and value != 0:
                counter += 1
                color = pallette.get(value, 'gray')
                # In case where j coupled to i, reverse the source - target
                if value == 4:
                    G.add_edge(target, source, color=color)
                G.add_edge(source, target, color=color)

    # Compute the degree of each node
    node_degrees = dict(G.degree())

    # Create a subgraph with nodes that have non-zero degree
    subgraph_nodes = [node for node, degree in node_degrees.items() if degree > 0]

    # Remove routes that include exchange reactions
    if model is not None and remove_exchange_routes:
        exchange_reactions = []
        for r in model.reactions:
            if 'e' in r.compartments or 'e0' in r.compartments:
                exchange_reactions.append(r.id)
        tmp = subgraph_nodes.copy()
        for node in tmp:
            if node in exchange_reactions:
                subgraph_nodes.remove(node)

    # Build a graph based on the subgraph nodes
    qfca_graph = G.subgraph(subgraph_nodes)

    if model:
        # graph's nodes are reactions, we ll add model name and stoichiometry as attributes
        for node in qfca_graph.nodes:
            r = model.reactions.get_by_id(node)
            reactants = [i.name for i in r.reactants]
            products = [i.name for i in r.products]
            qfca_graph.nodes[node]['rxn_name'] = r.name
            qfca_graph.nodes[node]['rxn_reactants'] = ';'.join(reactants)
            qfca_graph.nodes[node]['rxn_products'] = ';'.join(products)

    return qfca_graph


def samples_on_qfca(qfca_graph, samples):
    """
    ongoing
    """
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import KMeans
    from sklearn import metrics

    graph_rxns = list(qfca_graph.nodes)
    # df = np.zeros( [len(list(qfca_graph.nodes)), samples.shape[1] ])
    samples_with_graph_rnxs = samples.loc[graph_rxns]

    # Assuming df is your DataFrame containing continuous values

    # Initialize a list to store inertias
    inertias = []

    # Define the range of k values you want to test
    k_range = range(1, 11)  # You can adjust this range as needed

    # Calculate inertia for each k
    for k in k_range:
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(samples_with_graph_rnxs)
        inertias.append(kmeans.inertia_)


# %% Util functions: flux balance analysis
def producing_or_consuming_a_met(model, reaction_id, metabolite_id):
    """
    Returns whether a metabolite is being produced or consumed when model is optimized.

    Args:
        model (cobra.Model)
        reaction_id (str)
        metabolite_id (str)

    Returns:
        'producing | consuming' (str)
    """
    r = model.reactions.get_by_id(reaction_id)
    met = model.metabolites.get_by_id(metabolite_id)
    flux_value = r.summary().to_frame()["flux"].item()

    if met not in r.reactants and met not in r.products:
        return ValueError("the metabolite you are asking is not part of the reaction you gave.")

    if (flux_value > 0 and met in r.products) or ( flux_value < 0 and met in r.reactants):
        return "producing"
    elif (flux_value > 0 and met in r.reactants) or (flux_value < 0 and met in r.products):
        return "consuming"


def get_reactions_producing_a_met(model, metabolite_id):
    """
    Get reactions that produce a specific metabolite.

    Notes:
        Exloits the producing_or_consuming_a_met() function

    Args:
        model (cobra.Model)
        metabolite_id (str)

    Returns:
        list of reactions (cobra.Reaction)
    """
    rxns = []
    for reaction in model.metabolites.get_by_id(metabolite_id).reactions:
        if producing_or_consuming_a_met(model, reaction_id=reaction.id, metabolite_id=metabolite_id) == "producing":
            rxns.append(reaction)
    return rxns


def trace_path_with_backtracking_iterative(model, start_reaction_id, target_reaction_id, ignore_mets=None):
    """
    Trace a path from the start reaction to the target reaction through reactants,
    using an iterative DFS approach, with backtracking when exchange reactions are encountered.

    Args:
        model (cobra.Model): COBRApy model object.
        start_reaction_id (str): The ID of the starting reaction.
        target_reaction_id (str): The ID of the target reaction.
        ignore_mets: A list of metabolite IDs to ignore during tracing (optional).

    Returns:
        keep_rxns: A set of reactions needed to go from start to target, avoiding dead-end exchange reactions.
    """

    def is_exchange_reaction(reaction):
        """Check if a reaction is an exchange reaction."""
        return len(reaction.reactants) == 0 or len(reaction.products) == 0 or reaction.id.startswith("EX_")

    # Initialize structures for DFS
    visited_reactions = set()  # Tracks all reactions visited to prevent revisiting
    stack = [(model.reactions.get_by_id(start_reaction_id), [])]

    valid_paths = []  # List to store valid paths that do not dead-end at exchange reactions
    dead_end_paths = set()  # Tracks reactions leading to dead-end exchange reactions

    keep_rxns = set()  # Final set of reactions in valid paths

    while stack:
        current_reaction, path = stack.pop()

        # If we reached the target, store the path and continue exploring for other valid paths
        if current_reaction.id == target_reaction_id:
            valid_paths.append(path + [current_reaction.id])
            continue

        # Skip if the current reaction is an exchange reaction
        if is_exchange_reaction(current_reaction):
            dead_end_paths.add(current_reaction.id)
            continue

        visited_reactions.add(current_reaction.id)
        path.append(current_reaction.id)  # Track the path

        # Explore the reactants of the current reaction
        if current_reaction.summary().to_frame()["flux"].item() > 0:
            reactants = current_reaction.reactants
        else:
            reactants = current_reaction.products

        found_valid_branch = False  # Track if any valid branches are found from the current reaction
        for metabolite in reactants:

            if metabolite.id in ignore_mets:
                continue

            # Find reactions that produce this metabolite
            rxns = get_reactions_producing_a_met(model, metabolite_id=metabolite.id)
            for reaction in rxns:
                flux_value = reaction.summary().to_frame()["flux"].item()
                if flux_value == 0:
                    continue
                if reaction.id not in visited_reactions and reaction.id not in dead_end_paths:
                    # If we find a valid branch, continue exploring it
                    stack.append((reaction, path.copy()))  # Add new reaction to stack with the updated path
                    found_valid_branch = True
                    keep_rxns.add(reaction.id)

        # If no valid branch was found, mark this reaction as a dead-end (exchange reaction or no valid paths)
        if not found_valid_branch:
            dead_end_paths.add(current_reaction.id)

    # If no valid path found, return an empty set
    if not valid_paths:
        return set()

    # Return the reactions that are part of any valid path
    keep_rxns.add(start_reaction_id)
    return keep_rxns


