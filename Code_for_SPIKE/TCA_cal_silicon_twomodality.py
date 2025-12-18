import h5py
import tensortools as tt
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from joblib import Parallel, delayed
datafile = h5py.File('/mnt/Share/yuelp/TCAdecomposition/SPKdata.mat')
datamatrix = datafile['spike_normalized_laser_elec']
dataspike = np.array(datamatrix)
#dataspike = dataspike[10:, :, :]
rankrange = range(1, 31)
def TCAcal(rank,dataspike,outputfilename):
    ensemble = tt.Ensemble(fit_method='ncp_hals')
    ensemble.fit(dataspike, ranks=rank, replicates=10)
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
    with h5py.File(outputfilename, 'a') as file:
        file.create_dataset(str('factor_neuron_'+str(rank)), data=np.array(d))
        file.create_dataset(str('factor_trial_'+str(rank)), data=np.array(b))
        file.create_dataset(str('factor_time_' + str(rank)), data=np.array(c))
        file.create_dataset(str('factor_modality_' + str(rank)), data=np.array(a))
        file.create_dataset(str('factor_lam_' + str(rank)), data=np.array(lam))
        file.create_dataset(str('reconerror_' + str(rank)), data=np.array(reconerror))
        file.create_dataset(str('similarity_' + str(rank)), data=np.array(similarity))
        ...
    ...

Parallel(n_jobs=10)(delayed(TCAcal)(i,dataspike,'/mnt/Share/yuelp/TCAdecomposition/tensor_silicon_laser_elec.h5') for i in rankrange)


