import h5py
import tensortools as tt
import numpy as np
from joblib import Parallel, delayed
# tca calculate from the sampled datamatrix

def TCAcal(rank, dataspike, filename, iter):
    ensemble = tt.Ensemble(fit_method='ncp_hals')
    ensemble.fit(dataspike[iter, :, :, :], ranks=rank, replicates=10)
    a = list()
    b = list()
    c = list()
    d = list()
    lam = list()
    reconerror = list()
    similarity = list()
    for j in list(range(0, 10)):
        a.append(ensemble.results[rank][j].factors.factors[0])
        b.append(ensemble.results[rank][j].factors.factors[1])
        c.append(ensemble.results[rank][j].factors.factors[2])
        d.append(ensemble.results[rank][j].factors.factors[3])
        lam.append(ensemble.results[rank][j].factors.component_lams())
        reconerror.append(ensemble.results[rank][j].obj)
        similarity.append(ensemble.results[rank][j].similarity)
        ...
    ...
    with h5py.File(filename,'a') as file:
        file.create_dataset(str('factor_time_'+str(iter)), data=np.array(c))
        file.create_dataset(str('factor_trial_'+str(iter)), data=np.array(b))
        file.create_dataset(str('factor_neuron_' + str(iter)), data=np.array(d))
        file.create_dataset(str('factor_modality_' + str(iter)), data=np.array(a))
        file.create_dataset(str('factor_lam_' + str(iter)), data=np.array(lam))
        file.create_dataset(str('reconerror_' + str(iter)), data=np.array(reconerror))
        file.create_dataset(str('similarity_' + str(iter)), data=np.array(similarity))
    ...


datafile = h5py.File('/mnt/Share/yuelp/TCAdecomposition/SPKdata.mat')
datamatrix = datafile['spike_normalized_resample_laser_elec']
dataspike = np.array(datamatrix)
# run tensor
rank = 10 #
filename = '/mnt/Share/yuelp/TCAdecomposition/tensor_silicon_resample_laser_elec.h5'
Parallel(n_jobs=5)(delayed(TCAcal)(rank, dataspike, filename, iter) for iter in range(0, dataspike.shape[0]))

