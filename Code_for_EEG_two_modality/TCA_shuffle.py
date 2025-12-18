import tensortools as tt
import numpy as np
from joblib import Parallel, delayed
import sys
import h5py
import os


def TCA_shuffle(basepath, matfile, matvar, modality,outputfilename,rank):
    matfile = os.path.join(basepath, matfile)
    matvar = matvar
    outputfilename = os.path.join(basepath, outputfilename)
    modality = modality.split(',')
    rank = rank.split(',')
    datafile = h5py.File(matfile)
    datamatrix = np.array(datafile[matvar])  # time*frequency*channel*(trial*subject) time is np.linspace(-1,1.999,3000)
    #time = np.linspace(-1, 1.999, 3000)
    dataspike_slice = np.squeeze(datamatrix[:, list(np.array(modality, dtype='int') - 1), :, :, :, :])
    Parallel(n_jobs=30)(delayed(TCAcal)(rank, dataspike_slice, iter, outputfilename) for iter in range(0, datamatrix.shape[0]))


def TCAcal(rank, dataspike, iter, outputfilename):
    ensemble = tt.Ensemble(fit_method='ncp_hals')
    if np.ndim(dataspike)==6:
        ensemble.fit(np.squeeze(dataspike[iter, :, :, :, :, :]), ranks=int(rank[0]), replicates=10)
    else:
        ensemble.fit(np.squeeze(dataspike[iter, :, :, :, :]), ranks=int(rank[0]), replicates=10)
    ...
    a = list()
    b = list()
    c = list()
    d = list()
    e = list()
    lam = list()
    reconerror = list()
    similarity = list()
    for j in list(range(0, 10)):
        if np.ndim(dataspike) == 6:
            a.append(ensemble.results[rank][j].factors.factors[0])
            b.append(ensemble.results[rank][j].factors.factors[1])
            c.append(ensemble.results[rank][j].factors.factors[2])
            d.append(ensemble.results[rank][j].factors.factors[3])
            e.append(ensemble.results[rank][j].factors.factors[4])
        else:
            a.append(ensemble.results[rank][j].factors.factors[0])
            b.append(ensemble.results[rank][j].factors.factors[1])
            c.append(ensemble.results[rank][j].factors.factors[2])
            d.append(ensemble.results[rank][j].factors.factors[3])
        ...
        lam.append(ensemble.results[rank][j].factors.component_lams())
        reconerror.append(ensemble.results[rank][j].obj)
        similarity.append(ensemble.results[rank][j].similarity)
        ...
    ...
    with h5py.File(outputfilename, 'a') as file:
        if np.ndim(dataspike) == 6:
            try:
                file.create_dataset(str('factor_freq_' + str(rank)), data=np.array(c))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_trial_' + str(rank)), data=np.array(b))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_time_' + str(rank)), data=np.array(d))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_channel_' + str(rank)), data=np.array(e))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_modality_' + str(rank)), data=np.array(a))
            except:
                ...
            ...
        else:
            try:
                file.create_dataset(str('factor_freq_' + str(rank)), data=np.array(b))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_trial_' + str(rank)), data=np.array(a))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_time_' + str(rank)), data=np.array(c))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_channel_' + str(rank)), data=np.array(d))
            except:
                ...
            ...
        ...
        try:
            file.create_dataset(str('factor_lam_' + str(rank)), data=np.array(lam))
        except:
            ...
        ...
        try:
            file.create_dataset(str('reconerror_' + str(rank)), data=np.array(reconerror))
        except:
            ...
        ...
        try:
            file.create_dataset(str('similarity_' + str(rank)), data=np.array(similarity))
        except:
            ...
        ...
    ...


if __name__ == '__main__':
    rank = sys.argv[1]
    basepath = sys.argv[2]
    matfile = sys.argv[3]
    matvar = sys.argv[4]
    modality = sys.argv[5]
    outputfilename = sys.argv[6]
    TCA_shuffle(rank=rank, basepath=basepath, matfile=matfile, matvar=matvar, modality=modality, outputfilename=outputfilename)
...
