import h5py
import numpy as np

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

def gen_lattice_vectors(a,b,c,alpha,beta,gamma):
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
    latt_v =gen_lattice_vectors(*latt_list)
    yml_dict['lattice']['basis_vectors'] = ['{}'.format(v) for v in latt_v]
                                          
    yml_dict['excitations'] ={'type': 'DGSresolution'}
    UBi = np.linalg.inv(OL['UB'])
    u = np.dot(UBi,[0,0,1])
    v = np.dot(UBi,[0,1,0])
    yml_dict['orientation'] = {'u': u/np.abs(u).max(), 'v': v/np.abs(v).max() }
    yml_dict['shape'] = 'cylinder radius="12.5*mm" height="10.0*mm"' 
    yml_dict['temperature'] =  '0.3*K'
    
    

        



