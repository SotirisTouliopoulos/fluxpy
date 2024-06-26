
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
