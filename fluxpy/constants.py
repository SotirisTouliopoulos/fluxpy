import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap, to_rgba


def cmap():
    # Define your colors in hexadecimal format
    hex_colors = ['#9B32CD', '#DBD8DC', '#1FBDE7']

    # Convert hexadecimal colors to RGBA format
    rgba_colors = [to_rgba(color) for color in hex_colors]


    cmap_name = "my_cmap"
    n_bins = 100

    # sns.diverging_palette(220, 20, as_cmap=True, center='light', s=85, l=50)
    my_cmap = LinearSegmentedColormap.from_list(cmap_name, rgba_colors, N=n_bins)
    return my_cmap

TEST = "HELLO FRIEND"
CMAP = cmap()
