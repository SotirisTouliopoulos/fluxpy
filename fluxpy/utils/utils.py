""" fluxpy's main scope: provide functionalities to handle metabolic modeling analysis objects"""


import json
import cobra
import math
import numpy as np
import pandas as pd
import networkx as nx
from enum import Enum
from typing import List
from mergem import merge
from ..constants import *


# %% Inner functions
NestedDataFrameType = pd.DataFrame
all_namespaces = ["chebi", "metacyc", "kegg", "reactome", "metanetx", "hmdb", "biocyc", "bigg", "seed", "sabiork", "rhea"]

class _conversions(Enum):
    REACTION_NAMES = "reactionNames"
    REACTION_IDS = "reactionIds"
    METABOLITE_NAMES = "metaboliteNames"
    METABOLITE_IDS = "metaboliteIds"


def _convert_list_to_binary(lst):
    """
    Converts a ndarray to a binary decision based on the presence of at least a nonzero/True case.
    """
    if isinstance(lst, np.ndarray):
        if any(lst):
            return 1
        else:
            return 0
    else:
        return lst


def _convert_single_element_set(cell):
    if isinstance(cell, set) and len(cell) == 1:
        return next(iter(cell))  # Convert set to string
    return cell


# %% Util functions: parsing

def nonzero_fluxes(sol: pd.Series):
    """
    Returns nonzero fluxes from a pandas series
    """
    mask = sol.fluxes.to_numpy().nonzero()
    return sol.fluxes.iloc[mask]


def get_reactions_producing_met(model, met_id):
    """
    Returns a dataframe with the reaction ids of those consuming and producing met_id of a model.
    """

    solution = model.optimize()

    reactions_of_interest = model.metabolites.get_by_id(met_id).reactions

    reactions_producing_met = []
    reactions_consuming_met = []

    for reaction in reactions_of_interest:

        if solution[reaction.id] < 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_producing_met.append(reaction)
                    break

            for product in reaction.products:
                if met_id == product.id:
                    reactions_consuming_met.append(reaction)
                    break


        if solution[reaction.id] > 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_consuming_met.append(reaction)
                    break
            for product in reaction.products:
                if met_id == product.id:
                    reactions_producing_met.append(reaction)
                    break

    df = pd.DataFrame({'consuming_{met_id}': reactions_consuming_met, 'producing_{met_id}': reactions_producing_met})

    return df


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


def mapNamespaceWithModelSEEDDatabase(seed_compounds: List[str], path_to_modelSEEDDatabase: str, annotations_to_keep=['BiGG', 'BiGG1']):
    """
    This function makes use of the ModelSEEDpy and the ModelSEEDDatabase libraries to map a list of ModelSEED compounds
    to their annotations returning a pandas dataframe with those

    seed_compouds -- a list with seed compound ids
    path_to_modelSEEDDatabase -- path to the modelSEEDDatabase. To get it: git clone https://github.com/ModelSEED/ModelSEEDDatabase.git
    """
    from modelseedpy.biochem import from_local
    import gc  ## gc (garbage collection) module to force garbage collection.

    if path_to_modelSEEDDatabase is None:
        raise ValueError("path_to_modelSEEDdatabase is required.")
    try:
        modelseed_local = from_local(path_to_modelSEEDDatabase)
    except:
        raise ValueError("The modelSEEDDatabase provided was not correct. Please give correct path. For example: /Users/workspace/data/ModelSEEDDatabase/")
    dic = {}
    for medium_seed_met in seed_compounds:
        for met in modelseed_local.compounds:
            if met.id == medium_seed_met:
                found_annotations = {bigg: (bigg in met.annotation) for bigg in annotations_to_keep}
                if any(found_annotations.values()):
                    tmp = {annotation: met.annotation[annotation] if found else None for annotation, found in found_annotations.items()}
                    dic[medium_seed_met] = tmp
                else:
                    dic[medium_seed_met] = None
                break
    for key, value in dic.items():
        if value is None:
            dic[key] = {}
        for subkey in annotations_to_keep:
            dic[key].setdefault(subkey, None)
    t = pd.DataFrame.from_dict(dic, orient="index")
    df = t.applymap(_convert_single_element_set)
    del modelseed_local
    gc.collect()
    return df


def map2namespace(compounds, from_namespace="seed", to_namespace="bigg"):
    """
    compounds -- a list
    from_namesapce --
    to_namespace --
    """
    seed2metanex = pd.read_json(SEED2MNX)
    bigg2metanex = pd.read_json(BIGG2MNX)
    if from_namespace == "seed":
        seed_metanex = seed2metanex.loc[compounds].to_dict()["metanex_id"]
        metanex_ids = list(set(seed_metanex.values()))
        bigg_metanex = bigg2metanex[bigg2metanex["metanex_id"].isin(metanex_ids)]
        bigg_metanex = bigg_metanex["metanex_id"].to_dict()
        mapped_dict = {
            seed_key: next((bigg_key for bigg_key, bigg_value in bigg_metanex.items() if bigg_value == seed_value), None)
            for seed_key, seed_value in seed_metanex.items()
        }
        data_for_df = {'seed_metanex': list(mapped_dict.keys()), 'bigg_metanex': list(mapped_dict.values())}
    elif from_namespace == "bigg":
        bigg_metanex = bigg2metanex.loc[compounds].to_dict()["metanex_id"]
        metanex_ids = list(set(seed_metanex.values()))
        seed_metanex = seed2metanex[seed2metanex["metanex_id"].isin(metanex_ids)]
        seed_metanex = seed_metanex["metanex_id"].to_dict()
        mapped_dict = {
            bigg_key: next((seed_key for seed_key, seed_value in seed_metanex.items() if seed_value == bigg_value), None)
            for bigg_key, bigg_value in bigg_metanex.items()
        }
        data_for_df = {'bigg_metanex': list(mapped_dict.keys()), 'seed_metanex': list(mapped_dict.values())}
    df = pd.DataFrame(data_for_df)
    return df


# %% Util functions: flux analysis
def parse_qfca_output(qfca_output, model=None, remove_exchange_routes=True, exclude_biomass=True, format="csv"):
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

    # Add edges based on the values in the dataframe
    biomass_reaction = str(model.objective.expression).split()[0].split("*")[-1]
    for source in df.index:
        if exclude_biomass:
            if source == biomass_reaction:
                continue
        for target in df.index:
            if exclude_biomass:
                if target == biomass_reaction:
                    continue
            value = df.loc[source, target]
            if source != target and value != 0:  # Ignore cases where i = j or value is 0
                color = {1: 'black', 2: 'blue', 3: 'red', 4: 'green'}.get(value, 'gray')
                """
                1 - fully coupled reactions
                2 - partially coupled reactions
                3 - reaction i is directionally coupled to reaction j - red
                4 - reaction j is directionally coupled to reaction i - green
                """
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


# plt.plot(k_range, silhouette_scores, marker='o')
# plt.xlabel('Number of clusters (k)')
# plt.ylabel('Silhouette Score')
# plt.title('Silhouette Score for Optimal k')
# plt.show()

# %% Util classes

class Models:
    """
    Load a list of sbml models with cobra as attributes of a class.
    """
    def __init__(self, model_list: List[str], model_names: List[str] = None):
        if model_names is None:
            model_names =  ["".join(model.split("/")[-1].split(".")[:-1]) for model in model_list]
        for model_name, model_file in zip(model_names, model_list):
            setattr(self, model_name, cobra.io.read_sbml_model(model_file))


class NameIdConverter:
    """
    A class to switch indexes of a dataframe from ids to names and vice versa for metabolites and reactions of a model.

    model -- a cobra model
    df -- pandas dataframe to which the conversion will be applied
    convert -- current namespace of the df; values from ['reactionNames', reactionIds, 'metaboliteNames', 'metaboliteIds']
    to -- namespace to switch to; values as in `convert`
    """
    def __init__(self, model, df: pd.DataFrame, convert: str, to: str, init_only=False):

        try:
            convert, to = _conversions(convert), _conversions(to)
        except:
            raise ValueError(f"Both convert and to need to be among {[x.value for x in _conversions.mro()]}")

        met_map_df = pd.DataFrame([(met.id, met.name) for met in model.metabolites], columns=["id", "name"])
        react_map_df = pd.DataFrame([(react.id, react.name) for react in model.reactions], columns=["id", "name"])
        self.met_map_df = met_map_df
        self.react_map_df = react_map_df
        self.convert = convert
        self.to = to
        self.model = model
        self.df = df.to_frame().copy()
        self.isIndex = isinstance(df, pd.Index)
        if not init_only:
            self.run_convert()


    def _run_convert(self):

        convert_cases = {
            (self.convert.REACTION_NAMES.value,   self.to.REACTION_IDS.value): self.reactNames2reactIds,
            (self.convert.REACTION_IDS.value,     self.to.REACTION_NAMES.value): self.reactIds2reactNames,
            (self.convert.METABOLITE_NAMES.value, self.to.METABOLITE_IDS.value): self.metNames2metIds,
            (self.convert.METABOLITE_IDS.value,   self.to.METABOLITE_NAMES.value): self.metIds2metNames
        }
        # Retrieve the method associated with the specific combination of convert and to values
        convert_res = convert_cases.get((self.convert.value, self.to.value))
        if convert_res is not None:
            # Execute the corresponding function
            s = convert_res()
        else:
            # Handle the case when the combination doesn't exist
            raise ValueError("Invalid combination of convert and to values")

        self.df = s
        return self.df


    def _metIds2metNames(self):

        self.met_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.met_map_df['name'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
            self.df = pd.DataFrame(self.df.index, index = self.df["id"])
            self.df.columns = ["name"]
        return self.df

    def _metNames2metIds(self):

        self.met_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.met_map_df['id'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
        return self.df


    def _reactIds2reactNames(self):

        self.react_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.react_map_df['name'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["name"]
        return self.df


    def _reactNames2reactIds(self):

        self.react_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.react_map_df['id'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
        return self.df


class CompareModels:
    """
    Using mergem <https://mergem.readthedocs.io/> to compare a pair of models and returns a dataframe with
    reactions unique in model1 and model2 and those they share.

    Key arguments:\n
    model_list -- list of paths to model files to compare \n
    trans_to_db -- mergem allows to map metabolite and reaction ids to a series of namespaces.
    """

    def __init__(self, model_list: List[str], trans_to_db=None):

        self.number_of_models = len(model_list)
        self.model_names = ["".join(model.split("/")[-1].split(".")[:-1]) for model in model_list]
        self.models = Models(model_list)
        if trans_to_db is not None and trans_to_db not in all_namespaces:
            raise ValueError(f"namespace needs to be among {all_namespaces}")
        self.namespace = trans_to_db

        # Run mergem
        self.res = merge(model_list, set_objective='merge', trans_to_db=self.namespace)

        # Create presence-absence dataframes
        res = self.res.copy()
        for lst in res["met_sources"]:
            for i in range(self.number_of_models):
                if i not in res["met_sources"][lst]:
                    res["met_sources"][lst][i] = np.nan
        for lst in res["reac_sources"]:
            for i in range(self.number_of_models):
                if i not in res["reac_sources"][lst]:
                    res["reac_sources"][lst][i] = np.nan

        metabolites_df = pd.DataFrame.from_dict(res["met_sources"], orient="index")
        metabolites_df = metabolites_df.applymap(lambda x: np.where(pd.notnull(x), 1, 0))
        self.metabolites_df = metabolites_df.applymap(lambda x: _convert_list_to_binary(x))
        self.metabolites_df.columns = self.model_names

        reactions_df = pd.DataFrame.from_dict(self.res["reac_sources"], orient="index")
        reactions_df = reactions_df.applymap(lambda x: np.where(pd.notnull(x), 1, 0))
        self.reactions_df = reactions_df.applymap(lambda x: _convert_list_to_binary(x))
        self.reactions_df.columns = self.model_names


    def unique_metabolites(self, model_name):
        """
        Returns:
        unique_mets -- a pd.Index with the metabolites that are only present in model_name and not in any other of those in model_list
        """
        filtered_df = self.metabolites_df[self.metabolites_df[model_name] == 1]
        unique_mets = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return unique_mets

    def unique_reactions(self, model_name):
        """
        Key arguments:
        model_name -- the name (column name) of the model for which unique reactions will be found

        Returns:
        model_unique_reactions -- a pd.Index with the reactions that are only present in the model_name and not in any other of those in model_list
        """
        filtered_df = self.reactions_df[self.reactions_df[model_name] == 1]
        unique_reacts = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return unique_reacts


    def compare_model_pair(self, base_model=None, model_to_compare=None):
        """
        Compare pairwise models

        Key arguments:
        base_model -- name of the model as in dataframe, i.e. filename without the extension
        model_to_compare -- name of the model as in dataframe, i.e. filename without the extension

        Returns:
        r1 -- reactions only present in base model
        r2 -- reactions only present in model_to_compare
        m1 -- metabolites only present in base_model
        m2 -- metabolites only present in model_to_compare
        """

        if base_model is None and model_to_compare is None and self.number_of_models == 2:
            base_model, model_to_compare = self.model_names[0], self.model_names[1]

        reactions_m1 = self.reactions_df[self.reactions_df[base_model] == 1].index
        reactions_m2 = self.reactions_df[self.reactions_df[model_to_compare] == 1].index
        r1 = reactions_m1.difference(reactions_m2)
        r2 = reactions_m2.difference(reactions_m1)

        metabolites_m1 = self.reactions_df[self.metabolites_df[base_model] == 1].index
        metabolites_m2 = self.reactions_df[self.metabolites_df[model_to_compare] == 1].index
        m1 = metabolites_m1.difference(metabolites_m2)
        m2 = metabolites_m2.difference(metabolites_m1)
        return r1, r2, m1, m2


    def shared_metabolites(self, models_subset=None):
        """
        Key arguments:
        models_subset -- list of models to get shared metabolites and reactions

        Returns:
        shared_reactions -- a pd.Index with rearctions present among all models of models_subset
        """
        if models_subset is None:
            models_subset = self.model_names
        df = self.metabolites_df[models_subset]
        shared_metabolites = df[df.eq(1).all(axis=1)]

        return shared_metabolites.index


    def shared_reactions(self, models_subset=None):
        """
        Key arguments:
        models_subset -- list of models to get shared metabolites and reactions

        Returns:
        shared_reactions -- a pd.Index with rearctions present among all models of models_subset
        """
        if models_subset is None:
            models_subset = self.model_names
        df = self.reactions_df[models_subset]
        shared_reactions = df[df.eq(1).all(axis=1)]

        return shared_reactions.index





    # def compare_sublist(self, sublist):

    #     return 1



        # self.shared_metabolites =
        # self.shared_reactions =
        # self.unique_reactions_model1 =
        # self.unique_reactions_model2 =
        # self.unique_metabolites_model1 =
        # self.unique_metabolites_model2 =


