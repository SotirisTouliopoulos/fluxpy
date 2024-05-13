import numpy as np
import pandas as pd
import math
from enum import Enum

class MyEnum(Enum):
    REACTION_NAMES = "reactionNames"
    REACTION_IDS = "reactionIds"
    METABOLITE_NAMES = "metaboliteNames"
    METABOLITE_IDS = "metaboliteIds"

NestedDataFrameType = pd.DataFrame

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
    Assuming a single medium comound changes at a time return FVA findings
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


class convert_names_ids:
    """
    A class to switch indexes of a dataframe from ids to names and vice versa for metabolites and reactions of a model.

    model -- a cobra model
    df -- pandas dataframe to which the conversion will be applied
    convert -- current namespace of the df; values from ['reactionNames', reactionIds, 'metaboliteNames', 'metaboliteIds']
    to -- namespace to switch to; values as in `convert`
    """
    def __init__(self, model, df: pd.DataFrame, convert: str, to: str ):

        try:
            convert, to = MyEnum(convert), MyEnum(to)
        except:
            raise ValueError("Both convert and to need to be among 'reactionNames', reactionIds, 'metaboliteNames' and 'metaboliteIds'")

        met_map_df = pd.DataFrame([(met.id, met.name) for met in model.metabolites], columns=["id", "name"])
        react_map_df = pd.DataFrame([(react.id, react.name) for react in model.reactions], columns=["id", "name"])
        self.met_map_df = met_map_df
        self.react_map_df = react_map_df
        self.convert = convert
        self.to = to
        self.model = model
        self.df = df.to_frame().copy()
        self.isIndex = isinstance(df, pd.Index)


    def run_convert(self):

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


    def metIds2metNames(self):

        self.met_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.met_map_df['name'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
            self.df = pd.DataFrame(self.df.index, index = self.df["id"])
            self.df.columns = ["name"]
        return self.df

    def metNames2metIds(self):

        self.met_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.met_map_df['id'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
        return self.df


    def reactIds2reactNames(self):

        self.react_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.react_map_df['name'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["name"]
        return self.df


    def reactNames2reactIds(self):

        self.react_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.react_map_df['id'].to_dict()
        self.df.index = self.df.index.map(id_to_name_mapping)
        if self.isIndex:
            self.df.columns = ["id"]
        return self.df


class compare_models:
    """
    Use mergem <https://mergem.readthedocs.io/> to compare a pair of models and returns a dataframe with
    reactions unique in model1 and model2 and those they share
    """
    def __init__(self, model_list ):
        import mergem
        number_of_models = len(model_list)
        self.model_names = ["".join(model.split("/")[-1].split(".")[:-1]) for model in model_list]
        self.res = mergem.merge(model_list, set_objective='merge')
        res = self.res.copy()
        for lst in res["met_sources"]:
            for i in range(number_of_models):
                if i not in res["met_sources"][lst]:
                    res["met_sources"][lst].insert(i, np.nan)
        for lst in res["reac_sources"]:
            for i in range(number_of_models):
                if i not in res["reac_sources"][lst]:
                    res["reac_sources"][lst].insert(i, np.nan)
        self.metabolites_df = pd.DataFrame.from_dict(res["met_sources"], orient="index")
        self.metabolites_df = self.metabolites_df.applymap(lambda x: 1 if pd.notnull(x) else 0)
        self.metabolites_df.columns = self.model_names
        self.reactions_df = pd.DataFrame.from_dict(self.res["reac_sources"], orient="index")
        self.reactions_df = self.reactions_df.applymap(lambda x: 1 if pd.notnull(x) else 0)
        self.reactions_df.columns = self.model_names

    def unique_metabolites_for(self, model_name):
        filtered_df = self.metabolites_df[self.metabolites_df[model_name] == 1]
        indices = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return indices

    def unique_reactions_for(self, model_name):
        filtered_df = self.reactions_df[self.reactions_df[model_name] == 1]
        indices = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return indices


    def compare_model_pair(self, base_model, model_to_compare):
        """
        Returns:
        r1 -- reactions only present in base model
        r2 -- reactions only present in model_to_compare
        """
        reactions_m1 = self.reactions_df[self.reactions_df[base_model] == 1]
        reactions_m2 = self.reactions_df[self.reactions_df[model_to_compare] == 1]
        r1 = reactions_m1.index.difference(reactions_m2.index)
        r2 = reactions_m2.index.difference(reactions_m1.index)
        return r1, r2


    def shared_metabolites_among(self, models_subset=None):

        if models_subset is None:
            models_subset = self.model_names




    # def compare_sublist(self, sublist):

    #     return 1



        # self.shared_metabolites =
        # self.shared_reactions =
        # self.unique_reactions_model1 =
        # self.unique_reactions_model2 =
        # self.unique_metabolites_model1 =
        # self.unique_metabolites_model2 =


"""
import fluxpy
import cobra
input_model_list = ['gf_compl_default_bf/s_infantis_gf_default_bf_compl_medium.sbml', 'nn_gf_default_bf_mvl3A/s_infantis_nn_gf_mvl3_def_bf.xml', 'base/s_infantis_default_bf_no_gf.sbml']
f = fluxpy.utils.compare_models(input_model_list)
g = f.unique_metabolites_for("s_infantis_nn_gf_mvl3_def_bf")
h = fluxpy.utils.convert_names_ids(f.res["merged_model"], g, "metaboliteIds", "metaboliteNames")
h.run_convert()

m = cobra.io.read_sbml_model(input_model_list[0])
sol = m.optimize()

o = fluxpy.utils.convert_names_ids(m, sol.shadow_prices, "metaboliteIds", "metaboliteNames")
"""
