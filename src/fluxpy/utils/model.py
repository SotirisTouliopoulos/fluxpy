""" Utilities to parse, convert and edit cobra models """

import cobra
import numpy as np
import pandas as pd
from typing import List
from enum import Enum
from typing import List
from mergem import merge
from ..constants import *
from .utils import _convert_list_to_binary, _convert_single_element_set
from cobra.io import load_json_model
from collections import defaultdict

# %% Inner functions
all_namespaces = ["chebi", "metacyc", "kegg", "reactome", "metanetx", "hmdb", "biocyc", "bigg", "seed", "sabiork", "rhea"]

class _Conversions(Enum):
    REACTION_NAMES = "reactionNames"
    REACTION_IDS = "reactionIds"
    METABOLITE_NAMES = "metaboliteNames"
    METABOLITE_IDS = "metaboliteIds"




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

    """
    def __init__(self, model, to_convert: pd.DataFrame, convert: str, to: str, init_only=False):
        """
        model -- a cobra model
        to_convert -- pandas dataframe to which the conversion will be applied
        convert -- current namespace of the df; values from ['reactionNames', reactionIds, 'metaboliteNames', 'metaboliteIds']
        to -- namespace to switch to; values as in `convert`
        """

        try:
            convert, to = _Conversions(convert), _Conversions(to)
        except:
            raise ValueError(f"Both convert and to need to be among {[x.value for x in _Conversions.mro()]}")

        met_map_df = pd.DataFrame([(met.id, met.name) for met in model.metabolites], columns=["id", "name"])
        react_map_df = pd.DataFrame([(react.id, react.name) for react in model.reactions], columns=["id", "name"])
        self.met_map_df = met_map_df
        self.react_map_df = react_map_df
        self.convert = convert
        self.to = to
        self.model = model
        if isinstance(to_convert, pd.DataFrame):
            self.to_convert = to_convert.copy()
        elif isinstance(to_convert, cobra.core.solution.Solution):
            self.to_convert = to_convert.to_frame().copy()
        else:
            df = pd.DataFrame(to_convert)
            self.to_convert = df.copy()
        self.isIndex = isinstance(to_convert, pd.Index)
        if not init_only:
            self._run_convert()


    def _run_convert(self):

        convert_cases = {
            (self.convert.REACTION_NAMES.value,   self.to.REACTION_IDS.value): self._react_names_to_react_ids,
            (self.convert.REACTION_IDS.value,     self.to.REACTION_NAMES.value): self._react_ids_to_react_names,
            (self.convert.METABOLITE_NAMES.value, self.to.METABOLITE_IDS.value): self._met_names_to_met_ids,
            (self.convert.METABOLITE_IDS.value,   self.to.METABOLITE_NAMES.value): self._met_ids_to_met_names
        }
        # Retrieve the method associated with the specific combination of convert and to values
        convert_res = convert_cases.get((self.convert.value, self.to.value))
        if convert_res is not None:
            # Execute the corresponding function
            s = convert_res()
        else:
            # Handle the case when the combination doesn't exist
            raise ValueError("Invalid combination of convert and to values")

        self.to_convert = s
        return self.to_convert


    def _met_ids_to_met_names(self):

        self.met_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.met_map_df['name'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
            self.to_convert = pd.DataFrame(self.to_convert.index, index = self.to_convert["id"])
            self.to_convert.columns = ["name"]
        return self.to_convert

    def _met_names_to_met_ids(self):

        self.met_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.met_map_df['id'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
        return self.to_convert


    def _react_ids_to_react_names(self):

        self.react_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.react_map_df['name'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["name"]
        return self.to_convert


    def _react_names_to_react_ids(self):

        self.react_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.react_map_df['id'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
        return self.to_convert


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




def map_namespace_to_ModelSEEDDatabase(seed_compounds: List[str], path_to_modelSEEDDatabase: str, annotations_to_keep=['BiGG', 'BiGG1']):
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


def map_to_namespace(compounds, from_namespace="seed", to_namespace="bigg"):
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


def _perform_recursive(in_met, reaction, bigg_model, visited_metabolites, BIGG_COFACTORS, inner_reactions):
    # Check if metabolite is not in BIGG_COFACTORS and not visited
    if (
        in_met.id not in BIGG_COFACTORS and
        in_met.id not in visited_metabolites
    ):
        print(">>", in_met.id, "from", reaction.build_reaction_string())
        # Recursively find inner reactions for this metabolite
        find_inner_reactions(in_met, bigg_model,
                            visited_metabolites, BIGG_COFACTORS, inner_reactions)

def find_inner_reactions(metabolite, bigg_model, visited_metabolites, BIGG_COFACTORS, inner_reactions):

    # Mark metabolite as visited
    visited_metabolites.add(metabolite.id)

    reactions_with_the_met = bigg_model.metabolites.get_by_id(metabolite.id).reactions
    has_irreversible_reaction = any(not reaction.reversibility for reaction in reactions_with_the_met)

    # Check reactions involving this metabolite
    for reaction in reactions_with_the_met:

        # Add reaction to inner_reactions if it's not already there
        if reaction not in inner_reactions:
            inner_reactions.add(reaction)

            # Traverse through metabolites in the reaction
            if not reaction.reversibility:
                if metabolite in reaction.reactants:
                    for in_metabolite in reaction.products:
                        _perform_recursive(in_metabolite, reaction, bigg_model, visited_metabolites, BIGG_COFACTORS, inner_reactions)

                else:
                    print("metabolite:", metabolite.id, " not as reactant in:", reaction.build_reaction_string())
            else:
                print("metabolite:", metabolite.id, "is involved in the reversible reaction:", reaction.build_reaction_string())
                if len(reactions_with_the_met) == 1 or has_irreversible_reaction is False:
                    if metabolite in reaction.reactants:
                        for in_metabolite in reaction.products:
                            _perform_recursive(in_metabolite, reaction, bigg_model, visited_metabolites, BIGG_COFACTORS, inner_reactions)
                    else:
                        for in_metabolite in reaction.reactants:
                            _perform_recursive(in_metabolite, reaction, bigg_model, visited_metabolites, BIGG_COFACTORS, inner_reactions)




bigg_cofactors = ['atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'udp[c]', 'udp_c', 'ump[c]', 'ump_c',
                  'amp_c', 'amp[c]',
                  'gdp[c]', 'gdp_c', 'gtp[c]', 'gtp_c',
                  'accoa_c', 'accoa[c]', 'coa_c', 'coa[c]',
                  'q8[c]', 'q8_c', 'q8h2_c', 'q8h2[c]', 'mqn8_c', 'mqn8[c]', 'mql8_c', 'mql8[c]', 'q8h2_c', 'q8h2[c]',
                  'actp[c]', 'actp_c',
                  'h2o_c', 'h2o[c]', 'h2o_e', 'h2o[e]',
                  'pi_e', 'pi[e]', 'pi_c', 'pi[c]', 'ppi[c]', 'ppi_c',
                  'pep_c', 'pep[c]',
                  'h_c', 'h[c]', 'h_e', 'h[e]',
                  'o2_c', 'o2[c]', 'o2_e', 'o2[e]',
                  'co2_c', 'co2[c]', 'co2_e', 'co2[e]',
                  'nadp_c', 'nadp[c]', 'nadph_c', 'nadph[c]', 'nad_c', 'nad[c]', 'nadh_c', 'nadh[c]',
                  'nadp_e', 'nadp[e]', 'nadph_e', 'nadph[c]', 'nad_e', 'nad[e]', 'nadh_e', 'nadh[e]',
                  'fadh2_c', 'fadh2[c]', 'fad_c', 'fad[c]',
                  'nh4_c', 'nh4[c]', 'nh4_e', 'nh4[e]',
                  'pyr[c]', 'pyr_c'
                ]

bigg_building_blocks = ['ala_L[c]', 'asp_L[c]', ' gln_L[c]', 'glu_L[c]', 'glu_L[c]', 
                        'ser_L[c]', 'trp_L[c]', 'met_L[c]', 'lys_L[c]', 'cyst_L[c]'
                       ]


def find_arrow_index():
    # load cobra model
    cobra_model = load_json_model("../../../ext_data/")
    # list with reactions names
    reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    # define arrow types
    arrow_types = ["-->", "<=>", "<--"]
    # define list to store indices
    arrow_indices = []
    
    # parse every reaction
    for reaction in reactions_ids:
        reaction_stoichiometry = str(cobra_model.reactions.get_by_id(reaction))
        
        for arrow_type in arrow_types:
            # find index of arrow
            try:
                arrow_index = reaction_stoichiometry.split().index(arrow_type)
            except:
                arrow_index = None
            finally:
                if arrow_index != None:
                    arrow_indices.append(arrow_index)
                    
    return arrow_indices

                    
def find_reactants_products_cofactors():
    
    # get the arrow indices from the corresponding function
    arrow_indices = find_arrow_index()
    
    # load cobra model
    cobra_model = load_json_model("../../../ext_data/")
    
    # list with reactions names
    reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    # define variable found in the reaction stoichiometry str
    stoichiometry_elements = ["-->", "<=>", "<--", "+"]
    
    # define lists to store reactants, products and cofactors of all reactions
    reactants_all_reactions = []
    products_all_reactions = []
    cofactors_all_reactions = []
    
    for i in range(len(reactions_ids)):
        
        reaction = reactions_ids[i]
        arrow_index = arrow_indices[i]
        
        # define lists to store reactants, products and cofactors of one reaction
        reactants_single_reaction = []
        products_single_reaction = []
        cofactors_single_reaction = []
        
        # convert stoichiometry to str
        reaction_stoichiometry = str(cobra_model.reactions.get_by_id(reaction))

        # iterate through every element of the stoichiometry str, avoiding the first which is the reaction id
        for j in range(1, len(reaction_stoichiometry.split()) , 1):
            element = reaction_stoichiometry.split()[j]
            
            # try convert stoichiometries from str to float to find coefficients
            try:
                element = float(element)
            except:
                pass
            
            # find reactants
            if  element not in bigg_cofactors and element not in bigg_building_blocks and \
                element not in stoichiometry_elements and type(element) != float and j <= arrow_index:
                reactants_single_reaction.append(element)
                
            # find products
            elif  element not in bigg_cofactors and element not in bigg_building_blocks and \
                element not in stoichiometry_elements and type(element) != float and j >= arrow_index:
                products_single_reaction.append(element)

            # find cofactors
            elif element in bigg_cofactors or element in bigg_building_blocks:
                cofactors_single_reaction.append(element)
                     
        reactants_all_reactions.append(reactants_single_reaction)
        products_all_reactions.append(products_single_reaction)
        cofactors_all_reactions.append(cofactors_single_reaction)
        
    return reactions_ids, reactants_all_reactions, products_all_reactions, cofactors_all_reactions


def find_reactions_combinations():
    
    reactions_ids, reactants_all_reactions, products_all_reactions, cofactors_all_reactions = find_reactants_products_cofactors()
    
    # store possible combinations of reactions with same reactants/products and different cofactors
    cofactor_specificity = []

    for i in range(len(reactions_ids)):
        for j in range(len(reactions_ids)):
            
            # avoid comparison of the same reaction
            if reactions_ids[i] != reactions_ids[j]:

                # avoid comparison with empty lists
                if len(reactants_all_reactions[i]) > 0 and len(reactants_all_reactions[j]) > 0 and \
                   len(products_all_reactions[i]) > 0 and len(products_all_reactions[j]) > 0:
                
                    # boolean to check if pairwise reactants/products are the same
                    identical_reactants = (set(reactants_all_reactions[i]) == set(reactants_all_reactions[j]))
                    identical_products = (set(products_all_reactions[i]) == set(products_all_reactions[j]))
                    
                    if identical_reactants == True and identical_products == True:
                        # boolean to check if pairwise cofactors are the same
                        identical_cofactors = (set(cofactors_all_reactions[i]) == set(cofactors_all_reactions[j]))
                        
                        if identical_cofactors == False:
                            cofactor_specificity.append((reactions_ids[i], reactions_ids[j]))
                            
                            # for debugging only
                            #print("Reaction 1:", reactions_ids[i], "| Reaction 2:", reactions_ids[j])
                            #print("Reaction 1 Reactants:", reactants_all_reactions[i], "| Reaction 2 Reactants:", reactants_all_reactions[j])
                            #print("Reaction 1 Products:", products_all_reactions[i], "| Reaction 2 Products:", products_all_reactions[j])                        
                            #print("Reaction 1 Co-factors:", cofactors_all_reactions[i], "| Reaction 2 Co-factors:", cofactors_all_reactions[j])

    unique_combinations = {tuple(sorted(pair)) for pair in cofactor_specificity}
    #print("Unique combinations of reactions:" , list(unique_combinations))
    
    return unique_combinations


def merge_connected_combinations():

    unique_combinations = find_reactions_combinations()
    
    # Step 1: Create a graph using an adjacency list representation
    # We'll use a defaultdict where each key is a node, and its value is a set of neighboring nodes
    graph = defaultdict(set)

    # Build the adjacency list by iterating through each pair in the input
    for x, y in unique_combinations:
        # Add y as a neighbor of x
        graph[x].add(y)
        # Add x as a neighbor of y (because the graph is undirected)
        graph[y].add(x)

    # At this point, for an input of this format: [('A', 'B'), ('A', 'C'), ('B', 'C'), ('D', 'K')]
    # the `graph` looks like:
    # {'A': {'B', 'C'}, 'B': {'A', 'C'}, 'C': {'A', 'B'}, 'D': {'K'}, 'K': {'D'}}

    # Step 2: Define a function to find connected components using Depth-First Search (DFS)
    def find_connected_components(graph):
        # A set to keep track of visited nodes
        visited = set()
        # A list to store the connected components
        components = []

        # Helper function to perform a recursive DFS
        def dfs(node, component):
            # Mark the current node as visited
            visited.add(node)
            # Add the node to the current component
            component.append(node)
            # Visit all unvisited neighbors of the current node
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, component)

        # Iterate through all nodes in the graph
        for node in graph:
            # If the node has not been visited, it's the start of a new connected component
            if node not in visited:
                # Create a new component (list) to hold connected nodes
                component = []
                # Perform DFS starting from this node
                dfs(node, component)
                # Add the completed component to the list of components
                components.append(component)

        # Return the list of connected components
        return components

    # Step 3: Find connected components in the graph
    groups = find_connected_components(graph)

    # At this point, `groups` will contain:
    # [['A', 'B', 'C'], ['D', 'K']]

    # Print the result
    print(groups)
