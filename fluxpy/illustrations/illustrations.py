"""Functionalities to visualize findings of various types of metabolic modeling analysis"""

import math
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from dash import Dash, html
import dash_cytoscape as cyto
import plotly.graph_objs as go
from ..constants import *

# """
# Envelope analysis
# """
def plot_prod_env_3D(v1: pd.Series, v2: pd.Series, v3: pd.Series, width=800, height=600):
    """
    This function takes as arguments 3 columns of the cobra production_envelope() result
    to return a 3D scatter plot of those.

    Keyword arguments:
    v1 -- flux vector for x-axis
    v2 -- flux vector for y-axis
    v3 -- flux vector for z-axis which is used for the weight values

    Usage:
    from cobra.io import load_model
    model = load_model("textbook")
    from cobra.flux_analysis import production_envelope
    prod_env = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
    prod_env.head(3)
    carbon_source	flux_minimum	carbon_yield_minimum	mass_yield_minimum	flux_maximum	carbon_yield_maximum	mass_yield_maximum	EX_glc__D_e	EX_o2_e
    0	EX_glc__D_e	0.0	0.0	0.0	0.000000	0.000000	0.000000	-10.0	-60.000000
    1	EX_glc__D_e	0.0	0.0	0.0	1.578947	0.052632	0.051748	-10.0	-56.842105
    2	EX_glc__D_e	0.0	0.0	0.0	3.157895	0.105263	0.103496	-10.0	-53.684211

    x=prod_env['EX_o2_e']
    y=prod_env['EX_glc__D_e']
    z=prod_env['flux_maximum']
    """
    # Create a trace
    trace = go.Scatter3d(
        x=v1, y=v2, z=v3,
        mode='markers',
        marker=dict(
            size=16,
            color=v3,                # Set color to the z-axis value
            colorscale='emrld',   # Choose a colorscale
            opacity=0.8,
            colorbar=dict(title=f'{v3.name}')
        )
    )

    # Create layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title=v1.name),
            yaxis=dict(title=v2.name),
            zaxis=dict(title=v3.name)
        ),
        margin=dict(l=0, r=0, b=0, t=0),  # Adjust margin to remove unnecessary space
        width=width,
        height=height,
        # paper_bgcolor='rgba(0,0,0,0)',
        # plot_bgcolor='rgba(0,0,0,0)'
        # font=dict(
        #     family="Courier New, monospace",
        #     size=18,
        #     color="RebeccaPurple",
        #     variant="small-caps",
        # )
    )

    # Create a figure
    fig = go.Figure(data=[trace], layout=layout)


    # Show the plot
    return fig


def plot_ranging_medium_compounds(model, dictionary_with_plots, dpi=500):
    """
    Plots the behavior of medium compounds in a metabolic model across different concentrations.

    Args:
        model (object): The metabolic model object.
        dictionary_with_plots (dict): A dictionary where keys are compound IDs, and values are tuples
                                      containing lists of concentrations and corresponding flux values to plot.
        dpi (int, optional): Dots per inch (resolution) for the saved figure. Default is 500.

    Returns:
        A .png file
    """
    plt.rcParams.update({'font.size': 4})
    fig, axs = plt.subplots(math.ceil(len(dictionary_with_plots)/6), 6)
    counter = 0
    x = 0
    y = 0
    for compound, plot in dictionary_with_plots.items():
        compund_name = model.metabolites.get_by_id(compound[3:]).name
        x_axis = []
        for conc in plot[0]:
            conc = round(conc, 3)
            x_axis.append(conc)

        axs[x, y].plot(x_axis, plot[1])
        axs[x, y].set_title(compund_name)
        if counter == 0:
            ax1 = axs[x, y]
        else:
            axs[x, y].sharex(axs[0, 0])

        counter += 1
        y += 1
        if y == 6:
            x += 1
            y = 0

    fig.text(0.001, 0.5, 'FBA optimizing biomass flux', va='center', rotation='vertical', fontsize = 6)
    fig.text(0.3, 0.001, 'medium ranging the flux of one medium compound at a time', va='center', fontsize = 6)
    fig.tight_layout()
    fig.savefig("medium_compounds_behavior.png", bbox_inches='tight', dpi=dpi)


def plot_nutrients_gradient(gradient, nutrients=None, threshold=0.2, width_size=12, height_size=12, save_fig=False, path="gradient.png", dpi=500):
    """
    For each nutrient on a gradient it returns differences for the minimum and maximum cases
    and plots a heatmap for each where reaction that are affected with the increase of the upper bound
    at least threshold*max_difference are displayed.

    gradient -- output of the utils.get_nutrients_gradient()
    nutrients -- list with the reaction ids to be considered
    threshold -- per centage of the max difference that will be used to filter differences
    width_size -- width of each subplot
    height_size -- height of eash subplot
    """
    gradient = gradient.copy()
    if nutrients is None:
        nutrients = list(gradient.index)
    gradient = gradient.loc[nutrients]

    fig, axs = plt.subplots(len(nutrients), 2, figsize=(width_size, height_size*len(nutrients)))

    all_min_dfs = []
    all_max_dfs = []

    for nutrient_counter, (_, cell) in enumerate(gradient.iterrows()):

        nutrient = nutrients[nutrient_counter]

        # Calculate differences in minimum and maximum values between consecutive cells
        min_diffs = []
        max_diffs = []
        prev_df = None
        main_df = cell.to_frame().T
        for col in main_df.columns:
            df = main_df[col].iloc[0]
            if prev_df is not None:
                min_diff = df['minimum'] - prev_df['minimum']
                max_diff = df['maximum'] - prev_df['maximum']
                min_diffs.append(min_diff)
                max_diffs.append(max_diff)
            prev_df = df

        # Convert differences into DataFrame format
        min_diffs_df = pd.DataFrame(min_diffs, index = main_df.columns[1:], columns=df.index).T
        max_diffs_df = pd.DataFrame(max_diffs, index = main_df.columns[1:], columns=df.index).T

        all_min_dfs.append(min_diffs_df) ; all_max_dfs.append(max_diffs_df)

        # Calculate the maximum absolute difference across all rows and columns
        max_abs_diff_for_min = min_diffs_df.abs().values.flatten().max()
        max_abs_diff_for_max = max_diffs_df.abs().values.flatten().max()

        # Filter out rows that do not meet the threshold
        min_diffs_df_filtered = min_diffs_df[(min_diffs_df.abs() >= threshold * max_abs_diff_for_min).any(axis=1)]
        max_diffs_df_filtered = max_diffs_df[(max_diffs_df.abs() >= threshold * max_abs_diff_for_max).any(axis=1)]

        # Plot heatmaps
        # MINIMUM
        is_zero_df = min_diffs_df_filtered.any()
        """
        Returns a df like:
        ub_50     False
        ub_100    False
        dtype: bool
        """
        if is_zero_df.any():
            try:
                ax = axs[nutrient_counter, 0]
            except:
                ax = axs[0]
            maxd = max(abs(min_diffs_df_filtered.values.min()), abs(min_diffs_df_filtered.values.max()))
            norm = mcolors.TwoSlopeNorm(vmin=-maxd, vcenter=0, vmax=maxd)
            sns.heatmap(min_diffs_df_filtered, ax=ax, cmap=CMAP, cbar=True, annot=True, norm=norm)
            ax.set_title(f'Gradient over {nutrient} - minimum values')
            ax.set_xlabel(f'Change of min flux compared to the previous ub of {nutrient}')
            ax.set_ylabel('Reactions')
        # MAXIMUM
        is_zero_df = max_diffs_df_filtered.any()
        if is_zero_df.any():
            try:
                ax = axs[nutrient_counter, 1]
            except:
                ax = axs[1]
            maxd = max(abs(max_diffs_df_filtered.values.min()), abs(max_diffs_df_filtered.values.max()))
            norm = mcolors.TwoSlopeNorm(vmin=-maxd, vcenter=0, vmax=maxd)
            sns.heatmap(max_diffs_df_filtered, ax=ax, cmap=CMAP, cbar=True, annot=True, norm=norm)
            ax.set_title(f'Gradient over {nutrient} - maximum values')
            ax.set_xlabel(f'Change of max flux compared to the previous ub of {nutrient}')
            ax.set_ylabel('Reactions')

    plt.tight_layout()

    if save_fig:
        plt.savefig(path, bbox_inches='tight', dpi=dpi)

    plt.show()


    return all_min_dfs, all_max_dfs


# give a list of metabolites
# give their production/consumption based on experiments
# give thei sampling means, fba, fva something value
# compare in silico predictions against experiments


# """
# Coupling Analysis
# """
def qcfa_subgraphs(H):
    """
    Creates and runs a Dash application to visualize the QCFA subgraphs using a network graph.

    Args:
        H (networkx.Graph): A NetworkX graph where nodes represent reactions and edges represent relationships between them.
                            The edge colors denote different types of couplings.

    Returns:
        A Dash application visualizing the QFCA-oriented graph.
    """
    # Create Dash app
    app = Dash(__name__)

    # Draw the network with edge colors
    pos = nx.random_layout(H, seed=666)  # positions for all nodes

    # Define the layout of the app
    nodes = [{'data': {'id': node, 'label': node}, 'position': {'x': pos[node][0], 'y': pos[node][1]}} for node in H.nodes()]
    edge_id_counter=0
    edges = [{'data': {'id': f"edge_{edge_id_counter + i}", 'source': u, 'target': v}, 'classes': H[u][v]['color']} for i, (u, v) in enumerate(H.edges())]
    elements = nodes + edges

    directed_edges_selector = ""
    for edge in edges:
        u = edge['data']['source']
        v = edge['data']['target']
        if H[u][v]['color'] != "black":
            edge['classes'] += " triangle"
            directed_edges_selector += "#"+edge['data']['id']+', '

    directed_edges_selector = directed_edges_selector[:-2]

    legend_layout = html.Div([
        html.H3("Edge Legend"),
        html.Ul([
            html.Li([
                html.Span("Grey: ", style={'color': 'grey'}), "Fully coupled reactions"
            ]),
            html.Li([
                html.Span("Light Sea Green: ", style={'color': 'LightSeaGreen'}), "Partially coupled reactions"
            ]),
            html.Li([
                html.Span("Dark Salmon: ", style={'color': 'DarkSalmon'}), "Reaction i is directionally coupled to reaction j"
            ]),
            html.Li([
                html.Span("Steel Blue: ", style={'color': 'SteelBlue'}), "Reaction j is directionally coupled to reaction i"
            ])

        ])
    ])


    app.layout = html.Div([
        legend_layout,
        cyto.Cytoscape(
            id='network-graph',
            style={'width': '100%', 'height': '100vh'},
            elements=elements,
            stylesheet=[
                    {
                        'selector': 'node',  # Apply style to nodes
                        'style': {
                            'content': 'data(label)'  # Display node labels
                        }
                    },
                    {
                        'selector': '.red',
                        'style': {
                            'background-color': 'DarkSalmon',
                            'line-color': 'DarkSalmon',
                            'target-arrow-shape': 'none'
                        }
                    },
                    {
                        'selector': '.green',
                        'style': {
                            'background-color': 'LightSeaGreen',
                            'line-color': 'LightSeaGreen'
                        }
                    },
                    {
                        'selector': '.blue',
                        'style': {
                            'background-color': 'SteelBlue',
                            'line-color': 'SteelBlue'
                        }
                    },

                    {
                        'selector': '.black',
                        'style': {
                            'background-color': 'grey',
                            'line-color': 'grey'
                        }
                    },

                    {
                        'selector': directed_edges_selector,
                        'style': {
                            'background-color': 'data(background-color)',
                            'line-color': 'data(background-color)',  # Set line color based on data
                            'target-arrow-color': 'data(background-color)',  # Set arrow color based on data
                            'line-style': 'solid',  # Set line style to solid
                            'curve-style': 'bezier',  # Set curve style to bezier
                            'target-arrow-shape': 'triangle'  # Set arrow shape to triangle
                        }
                    },
            ],
            layout={
                'name': 'cose'
            },
        )
    ])

    # Run the app
    app.run_server(debug=True)


"""
Sampling analysis
"""
def flux_cone_example(facets=5, radius=5, cone_height=5, width=1000, height=1000):
    """
    Generates a 3D plot of a flux cone with random points inside it.

    Args:
        facets (int): Number of facets (sides) for the cone surface.
        radius (float): Radius of the base of the cone.
        cone_height (float): Height of the cone.
        width (int): Width of the plot in pixels.
        height (int): Height of the plot in pixels.

    Returns:
        plotly.graph_objects.Figure: A Plotly figure object containing the 3D plot of the cone and the points inside it.
    """

    # Generate random points inside the cone
    n_points = 1000
    x_points, y_points, z_points = _generate_points_inside_cone(n_points, radius, cone_height)

    # Create the figure
    fig = go.Figure()

    # Add the cone surface with transparency
    for i in range(facets):
        theta = i * 2 * np.pi / facets
        next_theta = (i + 1) * 2 * np.pi / facets
        x = [0, radius * np.cos(theta), radius * np.cos(next_theta)]
        y = [0, radius * np.sin(theta), radius * np.sin(next_theta)]
        z = [0, cone_height, cone_height]
        fig.add_trace(go.Mesh3d(x=x, y=y, z=z, color='rgba(0, 0, 255, 0.1)'))  # Adjust transparency here

    # Add random points
    fig.add_trace(go.Scatter3d(x=x_points, y=y_points, z=z_points, mode='markers',
                            marker=dict(size=3, color='red')))

    # Update layout
    fig.update_layout(title='Flux Cone with Random Points',
                    scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                    width=width , height=height)

    return fig


"""
Internal routines
"""
def _generate_points_inside_cone(n_points, radius, height):
    """
    Generates random points inside a cone.

    Args:
        n_points (int): Number of random points to generate.
        radius (float): Radius of the base of the cone.
        height (float): Height of the cone.

    Returns:
        tuple: Three numpy arrays containing the x, y, and z coordinates of the points inside the cone.
    """
    # Generate random points in cylindrical coordinates
    r = radius * np.sqrt(np.random.rand(n_points))
    theta = np.random.rand(n_points) * 2 * np.pi
    z = np.random.rand(n_points) * height

    # Convert cylindrical coordinates to Cartesian
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Filter out points outside the cone
    valid_indices = (x**2 + y**2) <= (radius * (z / height))**2
    x_inside = x[valid_indices]
    y_inside = y[valid_indices]
    z_inside = z[valid_indices]

    return x_inside, y_inside, z_inside
