""" Functionalities to handle metabolic modeling analysis objects"""

import numpy as np
import pandas as pd
import cobra
from typing import Union, Optional
from ..constants import *

NestedDataFrameType = pd.DataFrame


# %% Util functions: parsing

def nonzero_fluxes(fluxes:  Union[cobra.Solution, pd.Series]) -> None:
    """
    Returns nonzero fluxes from a pandas Series.

    This function takes a pandas Series containing flux values and returns a new Series with only the nonzero fluxes.

    Args:
        sol (pd.Series, mandatory): A pandas Series representing flux values.

    Returns:
        pd.Series: A pandas Series containing only the nonzero flux values from the input Series.
    """
    if not isinstance(fluxes, cobra.Solution) and not isinstance(fluxes, pd.Series):
        return ValueError("Provide a cobra.Solution or a pandas.Series element as fluxes.")
    if isinstance(fluxes, cobra.Solution):
        fluxes = fluxes.fluxes
    mask = fluxes.to_numpy().nonzero()
    return fluxes.iloc[mask]



def get_rxns_producing_consuming_met(met_id, model: cobra.Model = None, flux_vector: pd.Series = None):
    """
    Returns a DataFrame with the reaction IDs of those consuming and producing a given metabolite in a model.

    This function optimizes the given metabolic model to identify reactions involving a specified metabolite.
    It categorizes reactions into those producing and those consuming the metabolite and returns this information
    in a DataFrame.

    Args:
        model (cobra.Model, mandatory): The metabolic model to be analyzed.
        met_id (str, optional): The ID of the metabolite of interest.
        flux_vector (pd.Series, optional):

    Returns:
        pd.DataFrame: A DataFrame with two columns:
                    - `consuming_{met_id}`: Reactions consuming the specified metabolite in flux vector provided
                    - `producing_{met_id}`: Reactions producing the specified metabolite in flux vector provided

    Examples:

        >>> import cobra
        >>> model = cobra.io.read_sbml_model('e_coli_core.xml')
        >>> df = get_reactions_producing_met(model, 'met_id')
        >>> print(df)
            consuming_met_id producing_met_id
        0  reaction1          reaction2

    Notes:
        - The function assumes the model is properly optimized and that each reaction's flux can be accessed via the solution object.
        - Reactions with a flux less than 0 are considered to be consuming reactants and producing products in the reverse direction.
        - Reactions with a flux greater than 0 are considered to be producing reactants and consuming products in the forward direction.
        - This applies for the specific medium provided with the model
    """
    if model is None or not isinstance(model, cobra.Model):
        return TypeError("model needs to be a cobra.Model object.")

    if flux_vector is not None and not isinstance(flux_vector, pd.Series):
        return TypeError("flux_vector needs to be a pd.Series object.")

    if flux_vector is None:
        flux_vector = model.optimize()

    reactions_of_interest = model.metabolites.get_by_id(met_id).reactions

    reactions_producing_met = []
    reactions_consuming_met = []

    for reaction in reactions_of_interest:

        if flux_vector[reaction.id] < 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_producing_met.append(reaction)
                    break

            for product in reaction.products:
                if met_id == product.id:
                    reactions_consuming_met.append(reaction)
                    break

        if flux_vector[reaction.id] > 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_consuming_met.append(reaction)
                    break
            for product in reaction.products:
                if met_id == product.id:
                    reactions_producing_met.append(reaction)
                    break

    cons_df = pd.DataFrame({
        f'consuming {met_id}': reactions_consuming_met,
    })
    prod_df = pd.DataFrame({
        f'producing {met_id}': reactions_producing_met
    })
    return (cons_df, prod_df)


"""
Internal routines
"""
def _convert_list_to_binary(lst):
    """
    Converts an ndarray to a binary decision based on the presence of at least one nonzero/True value.

    This function checks if the input `lst` is a NumPy ndarray. If so, it evaluates whether there is
    at least one nonzero (or True) element in the array. If such an element exists, the function returns 1;
    otherwise, it returns 0. If `lst` is not a NumPy ndarray, it returns the original `lst` unchanged.

    Parameters:
    lst (any): The input that might be a NumPy ndarray.

    Returns:
    int or any: Returns 1 if `lst` is a NumPy ndarray with at least one nonzero/True value;
                returns 0 if `lst` is a NumPy ndarray with all elements zero/False;
                returns the original `lst` if it is not a NumPy ndarray.
    """
    if isinstance(lst, np.ndarray):
        if any(lst):
            return 1
        else:
            return 0
    else:
        return lst


def _convert_single_element_set(cell):
    """
    Converts a set with a single element to that element.

    This function checks if the input `cell` is a set containing exactly one element.
    If so, it converts the set to that single element. If `cell` is not a set with
    a single element, it returns the original `cell` unchanged.

    Parameters:
    cell (any): The input that might be a set with a single element.

    Returns:
    any: The single element if `cell` is a set with one element; otherwise, the original `cell`.
    """
    if isinstance(cell, set) and len(cell) == 1:
        return next(iter(cell))  # Convert set to string
    return cell

