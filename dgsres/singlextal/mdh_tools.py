import h5py
import numpy as np
import yaml
from mcvine.workflow import singlextal as sx
import warnings

class slice(object):
    def __init__(self,name):
        self.name = name
        self.grid = grid()
        self.res_2d_grid = grid()
        self.fitting = fitting()

class grid(object):
    def __init__(self):
        self.qaxis = None
        self.Eaxis = None

class fitting(object):
    def __init__(self):
        self.rounds = 3
        self.gaussian2d_threshold = 0.5
        self. alpha_bounds = (-np.pi/2, np.pi/2)


def sample_from_MDH(fl_name):
    """ get the sample info from the MDH file"""
    OL = {}
    with h5py.File(fl_name,'r') as fh:
        smpl_pth = '/MDHistoWorkspace/experiment0/sample'
        OL_pth = smpl_pth+"/oriented_lattice"
        kys = list(fh[OL_pth])
        OL['name'] = fh[smpl_pth].attrs['name'].decode()
        for ky in kys:
            if ky.find('unit_cell')>=0:
                OL[ky.split('cell_')[1]] = fh['{}/{}'.format(OL_pth,ky)][:][0]
            if ky.find('orientation_matrix')>=0:
                OL['UB'] = fh['{}/{}'.format(OL_pth,ky)][:]
    return OL

def angles_from_MDH(fl_name):
    """get the angle info from each experiment in an MDH"""
    with h5py.File(fl_name,'r') as fh:
        jq = list(fh["/MDHistoWorkspace"].keys())
        explist = [i for i in jq if i.find('experiment')>=0]
        angles = np.zeros(len(explist))
        for idx,exp in enumerate(explist):
            angles[idx] = fh["/MDHistoWorkspace/{}/logs/omega/value".format(exp)][:].mean()
    angles = np.sort(angles)
    dangles = angles[1:]-angles[:-1]
    if len(np.unique(dangles))>1:
        warnings.warn("Warning the angles are not equally spaced")
    return sx.axis(min=angles[0], max=angles[-1],step=dangles[0])    

def slice_from_MDH(fl_name,slice_name):
    """ get the slice info from an MDH"""
    sl = slice(slice_name)
    with h5py.File(fl_name,'r') as fh:
        projection = fh['MDHistoWorkspace/experiment0/logs/W_MATRIX/value'][:]
        projection = projection.reshape((3,3)) # need to check that this reshapes the matrix correctly.
        data_shp = np.array(fh['MDHistoWorkspace/data/signal'].shape)[::-1] # the shape of the data the first dimension is the last item in the tuple thus why reversing the array
        singledims = data_shp==1
        if singledims.sum()<2:
            raise RuntimeError('Must be a slice or a cut not a volume')
        Qdims = np.where(np.invert(singledims[:3]))[0]  # An array of booleans for which Q dimensions vary
        Q_perp_dims = np.where(singledims[:3])[0]  #An array of booleans for which Q dimensions are fixed.
        
        if singledims[-1]:
            # check if constant E cut # not completed yet.
            qs = {}
            for idx in range(Qdims):
                qs[idx] = fh['MDHistoWorkspace/data/D{}'.format(Qdims[idx])]
                
        
        else: # it is a constant Q cut
            hkl0 = np.zeros(3)
            for dimnum in range(len(Q_perp_dims)):
                qtmp = fh['MDHistoWorkspace/data/D{}'.format(Q_perp_dims[dimnum])][:].mean()
                hkl0 += projection[Q_perp_dims[dimnum],:]*qtmp     
            hkl_projection = projection[Qdims[0]]
            Etmp = fh['MDHistoWorkspace/data/D3'][:]
            Evals =(Etmp[1:]+Etmp[:-1])/2
            Eaxis = sx.axis(min=Evals.min(), max=Evals.max(), step=Evals[1]-Evals[0])
            qtmp = fh['MDHistoWorkspace/data/D{}'.format(Qdims[0])]
            qvals =(qtmp[1:]+qtmp[:-1])/2
            qaxis = sx.axis(min=qvals.min(), max=qvals.max(),step=qvals[1]-qvals[0])
            sl.hkl0 = hkl0
            sl.hkl_projection = hkl_projection
            sl.grid.qaxis = qaxis
            sl.grid.Eaxis = Eaxis
            return sl
            

def gen_lattice_vectors(a,b,c,alpha,beta,gamma):
    """a, b, c lattice parameters in Angstroms
    alpha, beta,gamma lattice angles in degrees
    returns a list of the basis vectors"""
    v1 = np.array([a,0,0])
    bprojx = b*np.cos(np.radians(gamma))
    bprojy = np.sqrt(b*b- bprojx*bprojx)
    v2 = np.array([bprojx,bprojy,0])
    cprojx = c*np.cos(np.radians(beta))
    cprojy = (b*c*np.cos(np.radians(alpha))-v2[0]*cprojx)/v2[1]
    cprojz = np.sqrt(c*c - cprojx**2-cprojy**2)
    v3 = np.array([cprojx, cprojy,cprojz])
    return [v1,v2,v3]

def prep_for_yml(OL):
    """ transform OL dictionary for use in Yaml"""
    yml_dict = {'name': OL['name'],'chemical_formula':OL['name']}
    latt_list = [OL['a'],OL['b'],OL['c'],OL['alpha'],OL['beta'],OL['gamma']]
    yml_dict['lattice'] = {'constants': '{},{},{},{},{},{}'.format(*latt_list) }
    latt_v = gen_lattice_vectors(*latt_list)
    yml_dict['lattice']['basis_vectors'] = [['{:0.4f}'.format(vi) for vi in vt] for vt in latt_v]                                    
    yml_dict['excitations'] ={'type': 'DGSresolution'}
    UBi = np.linalg.inv(OL['UB'])
    u = np.dot(UBi,[0,0,1])
    v = np.dot(UBi,[0,1,0])
    yml_dict['orientation'] = {'u': u/np.abs(u).max(), 'v': v/np.abs(v).max() }
    yml_dict['shape'] = 'cylinder radius="12.5*mm" height="10.0*mm"' 
    yml_dict['temperature'] =  '0.3*K'
    
    

        



