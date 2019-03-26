import numpy as np
from ..convolve2d import Grid

from ..convolve2d import Convolver as base
class Convolver(base):
    
    def __init__(
            self, grid, hkl_start, hkl_end, expansion_ratio, N_subpixels, res_func, res_range, transpose_res_matrix=True):
        base.__init__(self, grid, expansion_ratio, N_subpixels, res_func, res_range, transpose_res_matrix=transpose_res_matrix)
        self.hkl_start = hkl_start
        self.hkl_end = hkl_end
        qticks = grid.xaxis.ticks()
        newqticks = self.expanded_grid.xaxis.ticks()
        self.expanded_hkl_start, self.expanded_hkl_end \
            = expanded_hkl_range(
                hkl_start, hkl_end, qticks[0], qticks[-1], newqticks[0], newqticks[-1])
        return


def expanded_hkl_range(hkl_start, hkl_end, x_start, x_end, new_x_start, new_x_end):
    """the hkl start and end of the expanded range need to be calculated
    hkl_start and hkl_end corresponds to x_start and x_end
    now find the new_hkl_start and new_hkl_end that corresponds to new_x_start and new_x_end
    """
    vec = hkl_end-hkl_start
    univec = vec/(x_end-x_start)
    new_hkl_start = hkl_start + (new_x_start-x_start)*univec
    new_hkl_end = hkl_end + (new_x_end-x_end)*univec
    return new_hkl_start, new_hkl_end
