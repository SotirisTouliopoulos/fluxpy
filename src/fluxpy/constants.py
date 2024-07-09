"""Global variables for the fluxpy library"""


import seaborn as sns
import os
from matplotlib.colors import LinearSegmentedColormap, to_rgba


def _cmap():
    """
    Create a custom colormap.

    This function defines a custom colormap using a list of hexadecimal color codes,
    converts them to RGBA format, and then creates a LinearSegmentedColormap with
    100 bins.

    Returns:
        LinearSegmentedColormap: A custom colormap with the specified colors.
    """
    # Define your colors in hexadecimal format
    hex_colors = ['#9B32CD', '#DBD8DC', '#1FBDE7']

    # Convert hexadecimal colors to RGBA format
    rgba_colors = [to_rgba(color) for color in hex_colors]


    cmap_name = "my_cmap"
    n_bins = 100

    # sns.diverging_palette(220, 20, as_cmap=True, center='light', s=85, l=50)
    my_cmap = LinearSegmentedColormap.from_list(cmap_name, rgba_colors, N=n_bins)
    return my_cmap


CMAP = _cmap()

current_file_path = os.path.abspath(__file__)
parent_folder = os.path.dirname(current_file_path)
root_folder = os.path.dirname(parent_folder)
EXT_DATA = os.path.join(root_folder, "ext_data")
SEED2MNX = os.path.join(EXT_DATA, "seed2mnx.json")
BIGG2MNX = os.path.join(EXT_DATA, "bigg2mnx.json")

BIGG_COFACTORS = ['atp_c0', 'atp_c', 'adp_c', 'adp_c0',
                  'atp_c0', 'atp_c', 'adp_c', 'adp_c0',
                  'udp_c0', 'udp_c', 'ump_c0', 'ump_c',
                  'amp_c', 'amp_c0',
                  'gdp_c0', 'gdp_c', 'gtp_c0', 'gtp_c',
                  'accoa_c', 'accoa_c0', 'coa_c', 'coa_c0',  # acetyl-CoA
                  'q8_c0', 'q8_c', 'q8h2_c', 'q8h2_c0', 'mqn8_c', 'mqn8_c0', 'mql8_c', 'mql8_c0', 'q8h2_c', 'q8h2_c0',
                  'actp_c0', 'actp_c',
                  'h2o_c', 'h2o_c0', 'h2o_e', 'h2o[e]',
                  'pi_e', 'pi[e]', 'pi_c', 'pi_c0', 'ppi_c0', 'ppi_c',
                  'pep_c', 'pep_c0',
                  'h_c', 'h_c0', 'h_e', 'h[e]',
                  'o2_c', 'o2_c0', 'o2_e', 'o2[e]',
                  'co2_c', 'co2_c0', 'co2_e', 'co2[e]',
                  'nadp_c', 'nadp_c0', 'nadph_c', 'nadph_c0', 'nad_c', 'nad_c0', 'nadh_c', 'nadh_c0',
                  'nadp_e', 'nadp[e]', 'nadph_e', 'nadph_c0', 'nad_e', 'nad[e]', 'nadh_e', 'nadh[e]',
                  'fadh2_c', 'fadh2_c0', 'fad_c', 'fad_c0',
                  'nh4_c', 'nh4_c0', 'nh4_e', 'nh4[e]',
                  'pyr_c0', 'pyr_c'
                ]
BIGG_BUILDING_BLOCLS = ['ala_L_c0', 'asp_L_c0', ' gln_L_c0', 'glu_L_c0', 'glu_L_c0', 'ser_L_c0', 'trp_L_c0', 'met_L_c0', 'lys_L_c0', 'cyst_L_c0',

]

# Based on 10.1093/gigascience/giy021
MODELSEED_COFACTORS = [
    "cpd00004_c0",
    "cpd00005_c0",
    "cpd15353_c0",
    "cpd15499_c0",
    "cpd15561_c0",
    "cpd00097_c0",
    "cpd00982_c0",
    "cpd01270_c0",
    "cpd00002_c0",
    "cpd00038_c0",
    "cpd00052_c0",
    "cpd00062_c0",
    "cpd00068_c0",
    "cpd00115_c0",
    "cpd00241_c0",
    "cpd00356_c0",
    "cpd00357_c0",
    "cpd00358_c0",
    "cpd00530_c0",
    "cpd00977_c0",
    "cpd01775_c0"
]

QFCA_STYLE = os.path.join(root_folder, "fluxpy/data/style_qfca.json")

COMPLETE_MODEL = os.path.join(root_folder, "fluxpy/data/exchangeReactions.sbml")
