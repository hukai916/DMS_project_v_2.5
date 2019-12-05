"""
Selected packages from Enrich2 project;
mainly for plotting purposes.
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def recentered_cmap(cmap, vmin, vmax):
    """
    Rescale the diverging color map *cmap* such that the center color is at 0
    in the data. Returns the rescaled cmap.

    Based on http://stackoverflow.com/a/20528097

    *cmap* is the color map to re-center. Should be a diverging brewer
    palette.

    *vmin* is the minimum value in the data being plotted with the *cmap*
    (must be negative).

    *vmax* is the maximum value in the data being plotted with the *cmap*
    (must be positive).

    """
    # regular index to compute the colors
    reg_index = np.linspace(0.0, 1.0, 257)

    # shifted index to match the data
    centerpoint = 1 - vmax/(vmax + abs(vmin))
    #centerpoint = 0
    shift_index = np.hstack([
        np.linspace(0, centerpoint, 128, endpoint=False),
        np.linspace(centerpoint, 1, 129, endpoint=True)
    ])

    # re-map the colors at each index value using a dictionary
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    #print(reg_index)
    #print(shift_index)
    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        #print(ri,si)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    # save the dictionary as a color map
    newcmap = LinearSegmentedColormap("recentered", cdict)
    #print(cdict['red'])
    return newcmap
