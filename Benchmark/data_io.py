import numpy as np
import pandas as pd
import os

def read_and_process_AML(filepath):
    '''Read AML dataset from filepath, drop non-feature columns, save labels, arcsinh cofactor 5 transform the data
    in: string
    out: pandas dataframe'''
    datadf = pd.read_csv(filepath)
    datadf = datadf.loc[datadf['cell_type'] != 'NotDebrisSinglets']
    datadf = datadf.drop(axis='columns', labels=["Time", "Cell_length", "DNA1", "DNA2",
                                           "Viability", "file_number", "event_number", "subject"])
    datalabels = pd.DataFrame(datadf['cell_type'])
    datadf = datadf.drop(axis='columns', labels="cell_type")
    datadf = np.arcsinh((datadf-1)/5)
    return datadf, datalabels


def read_from_folder(filepath):
    labelfiles = os.listdir(filepath + '/Labels/')
    samplefiles = os.listdir(filepath + '/Samples/')
    labelfiles.sort()
    samplefiles.sort()
    samplesheap = pd.DataFrame()
    labelsheap = pd.DataFrame()

    for labelname, samplename in zip(labelfiles, samplefiles):
        samplesheap = samplesheap.append(pd.read_csv(filepath + '/Samples/' + samplename, header=None))
        labelsheap = labelsheap.append(pd.read_csv(filepath + '/Labels/' + labelname, header=None))
    return samplesheap, labelsheap

def read_and_process_BMMC(filepath):
    data = pd.read_csv(filepath)
    datalabels = pd.DataFrame(data['cell_type'])
    datalabels = datalabels[datalabels['cell_type'] != 'NotGated']
    datalabels = datalabels.reset_index().drop(axis='columns', labels='index')
    data = data[data['cell_type'] != "NotGated"]
    data = data.reset_index().drop(axis='columns', labels='index')
    data = data.drop(columns=['cell_type'])
    data = np.arcsinh((data-1)/5)
    return data, datalabels
