#%%Import libraries
import neo
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import os 
#import SpectralMUA_library
%matplotlib qt


#%%

import sys
sys.path.append(r"C:\Users\Javie\Documents\GitHub\SpectralEstimateOfMUA\Python code")
import SpectralMUA_library
#%%
#Path at laptop
LFPs_path = r"G:\My Drive\Manuscript SWO\Clean LFP exploration code\SMR\DLS"
#Path at desktop
#LFPs_path = r"D:\Javier\Papers\Dl-Dm Striatum\Data & Code\Dorsolateral data\SMRs"LFPs_list = [f for f in listdir(LFPs_path) if isfile(join(LFPs_path, f))]

LFPs_list = [f for f in listdir(LFPs_path) if isfile(join(LFPs_path, f))]

#%%
""" all_features = []
all_UDs = []
for rec in range(0,len(LFPs_list)):
    filename = os.path.join(LFPs_path+'/'+LFPs_list[rec])
    reader = neo.io.Spike2IO(filename=filename,try_signal_grouping=False)
    data = reader.read( lazy=False)[0]

    #data.segments[0].analogsignals[0].name #to know the channel name
    for i in range(0,len(data.segments[0].analogsignals)):
        if (data.segments[0].analogsignals[i].name == 'Vm'): Vm = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
        if (data.segments[0].analogsignals[i].name == 'M1'): motor = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
        if (data.segments[0].analogsignals[i].name == 'BFl'): somato = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
        if (data.segments[0].analogsignals[i].name == 'V1'): visual = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}

    #motor = np.asarray(data.segments[0].analogsignals[1])

    mua = SpectralMUA_library.computeMUA(motor['data'][:1939000],sampling_freq = motor['fs'])
    logMUA = np.log(mua)
    logMUA = SpectralMUA_library.smooth(logMUA,16) #80 ms smoothing as in Capone Cerebral Cortex 2020
    logMUA = logMUA[8:]
    parameters, Pfc_UDs_MUA = SpectralMUA_library.compute_UDs_logMUA(np.asarray(np.vstack((logMUA,logMUA))))
    #Note that the above line is extracting peaks from the MUA...
    #It can be correct, but I still need to validate it
    Pfc_UDs = []
    for i in range(0,len(Pfc_UDs_MUA)):
        Pfc_UDs.append(np.repeat(Pfc_UDs_MUA[i],5*motor['fs']/1000))
 """
#%%

MUAs = []
MUA_UDs = []
UDs = []
all_features = []
all_UDs = []
DLS = []
for rec in range(0,len(LFPs_list)):
    print('Trace nº '+str(rec))

    file = os.path.join(LFPs_path+'/'+LFPs_list[rec])
    temp_DLS = []
    temp_mua= []
    temp_mua_UD= []
    temp_UD= []
    temp_features = []
    reader = neo.io.Spike2IO(filename=file,try_signal_grouping=False)
    data = reader.read( lazy=False)[0]

    #data.segments[0].analogsignals[0].name #to know the channel name
    for i in range(0,len(data.segments[0].analogsignals)):
        if (data.segments[0].analogsignals[i].name == 'PF Ctex'): 
            FrA = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
            mua, mua_UD, UD, features = SpectralMUA_library.compute_logMUA_UDs(FrA['data'][:int(100*FrA['fs'])],FrA['fs'])
            temp_mua.append(mua)
            temp_mua_UD.append(mua_UD)
            temp_UD.append(UD)
            temp_features.append(features)
            temp_DLS.append(FrA['data'][:int(100*FrA['fs'])])
        if (data.segments[0].analogsignals[i].name == 'ORB'): 
            FrA = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
            mua, mua_UD, UD, features = SpectralMUA_library.compute_logMUA_UDs(FrA['data'][:int(100*FrA['fs'])],FrA['fs'])
            temp_mua.append(mua)
            temp_mua_UD.append(mua_UD)
            temp_UD.append(UD)
            temp_features.append(features)
            temp_DLS.append(FrA['data'][:int(100*FrA['fs'])])
            
        if (data.segments[0].analogsignals[i].name == 'M1'):
            motor = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
            mua, mua_UD, UD, features = SpectralMUA_library.compute_logMUA_UDs(motor['data'][:int(100*motor['fs'])],motor['fs'])
            temp_mua.append(mua)
            temp_mua_UD.append(mua_UD)
            temp_UD.append(UD)
            temp_features.append(features)
            temp_DLS.append(motor['data'][:int(100*motor['fs'])])

        if (data.segments[0].analogsignals[i].name == 'BFl'):
            somato = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
            mua, mua_UD, UD, features = SpectralMUA_library.compute_logMUA_UDs(somato['data'][:int(100*somato['fs'])],somato['fs'])
            temp_mua.append(mua)
            temp_mua_UD.append(mua_UD)
            temp_UD.append(UD)
            temp_features.append(features)
            temp_DLS.append(somato['data'][:int(100*somato['fs'])])

        if (data.segments[0].analogsignals[i].name == 'V1'): 
            visual = {'data': np.asarray(data.segments[0].analogsignals[i]), 'fs': np.asarray(data.segments[0].analogsignals[i].sampling_rate)}
            mua, mua_UD, UD, features = SpectralMUA_library.compute_logMUA_UDs(visual['data'][:int(100*visual['fs'])],visual['fs'])
            temp_mua.append(mua)
            temp_mua_UD.append(mua_UD)
            temp_UD.append(UD)
            temp_features.append(features)
            temp_DLS.append(visual['data'][:int(100*visual['fs'])])
    DLS.append(temp_DLS)
    MUAs.append(temp_mua)
    MUA_UDs.append(temp_mua_UD)
    all_UDs.append(temp_UD)
    all_features.append(temp_features)
    del temp_mua, temp_mua_UD, temp_UD, temp_DLS