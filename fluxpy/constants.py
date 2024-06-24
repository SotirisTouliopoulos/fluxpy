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

BIGG_COFACTORS = ['atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'udp[c]', 'udp_c', 'ump[c]', 'ump_c',
                  'amp_c', 'amp[c]',
                  'gdp[c]', 'gdp_c', 'gtp[c]', 'gtp_c',
                  'accoa_c', 'accoa[c]', 'coa_c', 'coa[c]',  # acetyl-CoA
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
BIGG_BUILDING_BLOCLS = ['ala_L[c]', 'asp_L[c]', ' gln_L[c]', 'glu_L[c]', 'glu_L[c]', 'ser_L[c]', 'trp_L[c]', 'met_L[c]', 'lys_L[c]', 'cyst_L[c]',

]


MODELSEED_COFACTORS = []


QFCA_STYLE = os.path.join(root_folder, "fluxpy/illustrations/style_qfca.json")