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

    res_func should be a function that returns a function.
    res_func(x, y) would return a new function f with which f(dxgrid, dygrid)
    would return the discretized resolution ellposid on the given grid.

    res_span: a tuple of x,y size of the PSF. For example, for a q,E slice, this would be
    (dq, dE) of the rectangle window for the resolution 
    """

    def __init__(self, grid, expansion_ratio, N_subpixels, res_func, res_range, transpose_res_matrix=True):
        self.grid = grid
        self.expansion_ratio = expansion_ratio
        self.N_subpixels = N_subpixels
        self.expanded_grid = self._expanded_grid()
        self.finer_expanded_grid = finer_grid(self.expanded_grid, N_subpixels)
        self.res_func = res_func
        dx, dy = res_range
        xstep, ystep = self.finer_expanded_grid.xaxis.step, self.finer_expanded_grid.yaxis.step
        Nx = int(np.ceil(dx/xstep)); hNx = int(Nx//2); Nx =hNx*2+1
        Ny = int(np.ceil(dy/ystep)); hNy = int(Ny//2); Ny =hNy*2+1
        dxs = np.arange(-hNx*xstep, (hNx+.5)*xstep, xstep)
        dys = np.arange(-hNy*ystep, (hNy+.5)*ystep, ystep)
        self.res_grids = np.meshgrid(dxs, dys)
        self.transpose_res_matrix = transpose_res_matrix
        return

    def convolve(self, img):
        """img is the data corresponding to the finer_expanded_grid
        """
        feg = self.finer_expanded_grid
        expected = (feg.xaxis.ticks().size, feg.yaxis.ticks().size)
        assert img.T.shape == expected, "image shape is %s, expecting %s" % (img.shape, expected)
        return convolve(img, self.psf)
    
    def psf(self, col, row):
        feg = self.finer_expanded_grid
        x = feg.xaxis.ticks()[col]
        y = feg.yaxis.ticks()[row]
        rt = self.res_func(x,y).ongrid(*self.res_grids)
        if self.transpose_res_matrix:
            rt = rt.T
        return rt

    def _expanded_grid(self):
        grid = self.grid
        ratio = self.expansion_ratio
        newxaxis = expand_axis_by_ratio(grid.xaxis, ratio)
        newyaxis = expand_axis_by_ratio(grid.yaxis, ratio)
        return Grid(newxaxis, newyaxis)

    
def convolve(image, psf, mode=''):
    """psf: a function that returns a PSF at any specific position in the image

    mode: just to be consistent with signatures of other convolve methods
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # for each pixel in image, replace it with the psf
    # be careful with the boundaries.
    # add everything together
    rows, cols = image.shape
    out = np.zeros(image.shape, dtype=float)
    for row in range(rows):
        # print row
        for col in range(cols):
            # total index
            i = row*cols + col
            if i%size != rank: continue
            # print "* cpu %s working on (%s, %s)" % (rank, col, row)
            # print row, col
            I = image[row, col]
            if I==0: continue
            spreaded = I*psf(col, row)
            psf_rows, psf_cols = spreaded.shape
            left = col - psf_cols//2
            right = left + psf_cols
            top = row - psf_rows//2
            bottom = top + psf_rows
            psf_left = psf_top = 0
            psf_right = psf_cols; psf_bottom = psf_rows
            # print top, bottom, left, right
            # print psf_top, psf_bottom, psf_left, psf_right
            if left < 0:
                shift = left
                left = 0
                psf_left -= shift
            if top < 0:
                shift = top
                top = 0
                psf_top -= shift
            if right > cols:
                shift = right - cols
                right = cols
                psf_right -= shift
            if bottom > rows:
                shift = bottom - rows
                bottom = rows
                psf_bottom -= shift
            # print top, bottom, left, right
            # print psf_top, psf_bottom, psf_left, psf_right
            out[top:bottom, left:right] += spreaded[psf_top:psf_bottom, psf_left:psf_right]
            continue
    tag = 100
    if rank != 0:
        comm.send(out, dest=0, tag=tag)
    else:
        for i in range(1, size):
            out1 = comm.recv(source=i, tag=tag)
            out+=out1
            continue
    out = comm.bcast(out, root=0)
    return out


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
