"""Routines to calculate stats on findings of metabolic modeling analysis"""


import numpy as np
from scipy.stats import linregress

def slope_in_prod_env(uptake_flux_vector, paused_uptake, objective_flux_vector, **kwargs):
    """
    Returns the slope of the line connecting the origin of coordinates and the
    point of the highest optimal value at the maximum of an uptake rate.
    For example, in case of growth rate and glucose,
    this slope corresponds to the biomass yield on glucose under anaerobic conditions,
    and in terms of linear programming, to the negative of the shadow price.

    Keyword arguments:
    uptake --
    objective --

    Returns:
    slope --  The slope of the line from the origin to uptake flux vector max
    shadow_price -- Shadow price of the corresponding uptake flux vector if a model provided
    """
    # Remove NaN values from x and y
    # valid_indices = ~np.isnan(uptake_flux_vector) & ~np.isnan(objective_flux_vector)
    # uptake_flux_vector_valid = uptake_flux_vector[valid_indices]
    # objective_flux_vector_valid = objective_flux_vector[valid_indices]

    # Calculate the slope using linear regression
    # slope, _, _, _, _ = linregress(uptake_flux_vector_valid, objective_flux_vector_valid)

    valid_indices_pair1 = ~np.isnan(uptake_flux_vector) & ~np.isnan(paused_uptake)
    valid_indices_pair2 = ~np.isnan(paused_uptake) & ~np.isnan(objective_flux_vector)
    valid_indices_pair3 = ~np.isnan(uptake_flux_vector) & ~np.isnan(objective_flux_vector)

    valid_indices = valid_indices_pair1 & valid_indices_pair2 & valid_indices_pair3

    uptake_flux_vector_valid = uptake_flux_vector[valid_indices]
    uptake_flux_vector_valid.reset_index(drop=True, inplace=True)
    paused_uptake_valid = paused_uptake[valid_indices]
    paused_uptake_valid.reset_index(drop=True, inplace=True)

    objective_flux_vector_valid = objective_flux_vector[valid_indices]
    objective_flux_vector_valid.reset_index(drop=True, inplace=True)

    # Find indices where z is 0
    indices = np.where(paused_uptake_valid == 0)[0]

    # Get the corresponding x values where z is 0
    x_values = uptake_flux_vector_valid[indices]

    # Find minimum and maximum x values where z is 0
    x_min = np.min(x_values)
    x_max = np.max(x_values)

    # Get corresponding z values for minimum and maximum x values
    z_min = np.min(objective_flux_vector_valid[np.where(uptake_flux_vector_valid == x_min)[0]])
    z_max = np.max(objective_flux_vector_valid[np.where(uptake_flux_vector_valid == x_max)[0]])

    print(x_min, z_min) ; print(x_max, z_max)

    slope = (z_max - z_min) / (x_max - x_min)

    print(f"Slope of the line for {uptake_flux_vector.name} (NaN values omitted):", slope)

    # Get shadow price
    model = kwargs.get('model', None)
    shadow_price = None
    if model is not None:
        sol = model.optimize()
        shadow_price = sol.shadow_prices.loc[uptake_flux_vector.name.replace("EX_", "")]
        print(f"Shadow price of the uptake reaction importing {uptake_flux_vector.name}:", shadow_price)

    return slope, shadow_price





