import h5py

def sample_from_MDH(fl_name):
    """ get the sample info from the MDH file"""
    with h5py.File(fl_name,'r') as fh:
        fh['/MDHistoWorkspace/experiment0/sample']
