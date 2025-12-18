import tensortools as tt
import numpy as np
from joblib import Parallel, delayed
import sys
import h5py
import os


def TCA_shuffle(basepath, matfile, matvar, outputfilename,rank):
    matfile = os.path.join(basepath, matfile)
    matvar = matvar
    outputfilename = os.path.join(basepath, outputfilename)
    datafile = h5py.File(matfile)
    datamatrix = datafile[matvar]  # time*frequency*channel*(trial*subject) time is np.linspace(-1,1.999,3000)
    time = np.linspace(-1, 1.999, 3000)
    datamatrix = datamatrix[:,:,:,:,np.logical_and(time >= -0.2, time <= 0.5)] # used slice to shuffle because long range time of shuffle is too large.
    datamatrix = np.array(datamatrix)
    Parallel(n_jobs=30)(delayed(TCAcal)(rank, datamatrix, iter, outputfilename) for iter in range(0, datamatrix.shape[0]))


def TCAcal(rank, dataspike, iter, outputfilename):
    ensemble = tt.Ensemble(fit_method='ncp_hals')
    ensemble.fit(dataspike[iter, :, :, :, :], ranks=rank, replicates=10)
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
        file.create_dataset(str('factor_freq_' + str(iter)), data=np.array(c))
        file.create_dataset(str('factor_trial_' + str(iter)), data=np.array(a))
        file.create_dataset(str('factor_time_' + str(iter)), data=np.array(d))
        file.create_dataset(str('factor_channel_' + str(iter)), data=np.array(b))
        file.create_dataset(str('factor_lam_' + str(iter)), data=np.array(lam))
        file.create_dataset(str('reconerror_' + str(iter)), data=np.array(reconerror))
        file.create_dataset(str('similarity_' + str(iter)), data=np.array(similarity))
    ...


if __name__ == '__main__':
    rank = sys.argv[1]
    basepath = sys.argv[2]
    matfile = sys.argv[3]
    matvar = sys.argv[4]
    outputfilename = sys.argv[5]
    TCA_shuffle(rank=int(rank), basepath=basepath, matfile=matfile, matvar=matvar, outputfilename=outputfilename)
...
