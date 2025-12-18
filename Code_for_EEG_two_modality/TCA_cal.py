import tensortools as tt
import numpy as np
from joblib import Parallel, delayed
import sys
import h5py
import os


def TCA_eval(basepath, matfile, matvar, modality, outputfilename):
    matfile = os.path.join(basepath, matfile)
    matvar = matvar
    modality = modality.split(',')
    outputfilename = os.path.join(basepath, outputfilename)
    datafile = h5py.File(matfile)
    datamatrix = datafile[matvar]  # time*frequency*channel*(trial*subject) time is np.linspace(-1,1.999,3000)
    dataspike = np.array(datamatrix)
    time = np.linspace(-1, 1.999, 3000)
    tindex = np.where(np.logical_and(time>=-0.2,time<=0.5))
    #dataspike_slice = np.squeeze(dataspike[list(np.array(modality,dtype='int')-1),:,:,:,:])
    for i in range(4):
        if not 'dataspike_slice' in locals():
            dataspike_slice = dataspike[:, i, list(np.array(modality, dtype='int')-1), :, :, :]
        else:
             dataspike_slice = np.concatenate((dataspike_slice, dataspike[:, i, list(np.array(modality,dtype='int')-1), :, :, :]), axis=0)
    ...
    if np.ndim(dataspike_slice) ==5:
        dataspike_slice = np.squeeze(dataspike_slice[:,:,:,:,tindex]) # trial*channel*freq*time
    else:
        dataspike_slice = np.squeeze(dataspike_slice[:,:,:,tindex]) # trial*modality*channel*freq*time
    ...
    rankrange = range(1, 31)
    Parallel(n_jobs=30)(delayed(TCAcal)(i, dataspike_slice, outputfilename) for i in rankrange)


def TCAcal(rank, dataspike, outputfilename):
    ensemble = tt.Ensemble(fit_method='ncp_hals')
    ensemble.fit(dataspike, ranks=rank, replicates=10)
    a = list()
    b = list()
    c = list()
    d = list()
    e = list()
    lam = list()
    reconerror = list()
    similarity = list()
    for j in list(range(0, 10)):
        if np.ndim(dataspike) == 5:
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
        if np.ndim(dataspike) == 5:
            try:
                file.create_dataset(str('factor_freq_' + str(rank)), data=np.array(d))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_trial_' + str(rank)), data=np.array(a))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_time_' + str(rank)), data=np.array(e))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_channel_' + str(rank)), data=np.array(c))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_modality_' + str(rank)), data=np.array(b))
            except:
                ...
            ...
        else:
            try:
                file.create_dataset(str('factor_freq_' + str(rank)), data=np.array(c))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_trial_' + str(rank)), data=np.array(a))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_time_' + str(rank)), data=np.array(d))
            except:
                ...
            ...
            try:
                file.create_dataset(str('factor_channel_' + str(rank)), data=np.array(b))
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
    basepath = sys.argv[1]
    matfile = sys.argv[2]
    matvar = sys.argv[3]
    modality = sys.argv[4]
    outputfilename = sys.argv[5]
    TCA_eval(basepath=basepath, matfile=matfile, matvar=matvar, modality=modality, outputfilename=outputfilename)
...
