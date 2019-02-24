import numpy as np

class Convolver:
    
    """convolve a 2D slice and a resolution fuction (position dependent)

    This class provides attributes and methods to make the convolution
    easier. One thing it needs to do is to consider adding edges.
    The grid given defines a window. For convolution, we will expand
    the window, then obtain data for the expanded window (which usually
    is a calculation result), and then perform the convolution,
    and at last trim it down to the requested size.

    N_subpixels: allow refining the slice grid. A good practice is to 
    perform calculation in a finer grid first before convert back to the
    original, coarser grid.
    """

    def __init__(self, grid, expansion_ratio, N_subpixels):
        self.grid = grid
        self.expansion_ratio = expansion_ratio
        self.N_subpixels = N_subpixels
        self.expanded_grid = self._expanded_grid()
        self.finer_expanded_grid = finer_grid(self.expanded_grid, N_subpixels)
        return

    def _expanded_grid(self):
        grid = self.grid
        ratio = self.expansion_ratio
        newxaxis = expand_axis_by_ratio(grid.xaxis, ratio)
        newyaxis = expand_axis_by_ratio(grid.yaxis, ratio)
        return Grid(newxaxis, newyaxis)

    
class Grid:

    def __init__(self, xaxis, yaxis):
        """xaxis and yaxis should be instances of dgsres.axis
        """
        self.xaxis = xaxis
        self.yaxis = yaxis
        return

    def __str__(self):
        return "Grid:\n - xaxis: %s\n - yaxis: %s" % (self.xaxis, self.yaxis)

    
def finer_grid(grid, N_subpixels):
    if isinstance(N_subpixels, tuple):
        Nx, Ny = N_subpixels
    else:
        Nx = Ny = N_subpixels
    x1 = finer_axis(grid.xaxis, Nx)
    y1 = finer_axis(grid.yaxis, Ny)
    return Grid(x1, y1)

def finer_axis(axis, N_subpixels):
    min = axis.min; max = axis.ticks()[-1] ; step = axis.step
    newstep = step/N_subpixels
    from . import axis
    return axis(min, max+newstep/2., newstep)

def expand_axis_by_ratio(axis, ratio):
    ticks = axis.ticks()
    min = ticks[0]; max = ticks[-1]; N = len(ticks); step = ticks[1]-ticks[0]
    n = int(np.ceil(N*ratio/2))
    newmin = min-n*step; newmax = max+n*step
    from . import axis
    return axis(newmin, newmax+step/2, step)


from . import axis
if '__str__' not in axis.__dict__:
    def axis_str(self):
        return 'axis in [%s, %s), step=%s' % (self.min, self.max, self.step)
    axis.__str__ = axis_str
