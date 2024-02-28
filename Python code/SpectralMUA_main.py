
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import loadmat
from os import listdir
from os.path import isfile, join
import os 
from scipy.signal import find_peaks 
from scipy.stats import zscore
%matplotlib qt
import SpectralMUA_library
#%% Functions to use


def ComputeD2USpeed_ind_UP(Transition):
      TimeWindow = np.arange(np.size(Transition))+1
      S,p = np.polyfit(TimeWindow,Transition,1)
      Slope = S
      return Slope*20000


def compute_parameters(dataset):
#      global dataset, parameters
    cell = []
    features = []
    UDs = []
    x = dataset
    for num in range(np.squeeze(np.shape(x[:,0]))):
        
            a = smooth(x[num,:], window_len=200)
            a = a[100:]
            a = np.transpose(a)
            UD = []
            D2U =[]
            U2D = []
            partition = int(np.round(np.size(a)/5))
            # for i in range(0,np.size(a)-partition,partition):
            for i in range(0,np.size(a),partition):
                temp = a[i:i+partition]
                thresh = np.mean(temp)+(0.5*np.std(temp))
                temp_UD = np.zeros(np.size(temp))
                for j in range(np.size(temp)):
                        if temp[j]>=thresh:
                            temp_UD[j] = 1
                UD.extend(temp_UD)
            UD = np.asarray(UD)
            for i in range(np.size(UD)-1):
                if (UD[i]==0 and UD[i+1]==1):
                    D2U.append(i)
                elif (UD[i]==1 and UD[i+1]==0):
                    U2D.append(i)
            D2U = np.asarray(D2U)
            U2D = np.asarray(U2D)
            for i in range(1,np.size(U2D)):
                temp = D2U[D2U>U2D[i]]
                if (np.size(temp)>0 and temp[0] - U2D[i]<250):
                    UD[U2D[i]-1:temp[0]+1] = 1
                
            D2U = []
            U2D = []
            for i in range(np.size(UD)-1):
                if (UD[i]==0 and UD[i+1]==1):
                    D2U.append(i)
                elif (UD[i]==1 and UD[i+1]==0):
                    U2D.append(i)
            D2U = np.asarray(D2U)
            U2D = np.asarray(U2D)
            for i in range(np.size(D2U)):
                temp = U2D[U2D>[D2U[i]]]
                if (np.size(temp)>0 and temp[0]-D2U[i]<20):
                    UD[D2U[i]:temp[0]+1] = 0
        
            D2U = []
            U2D = []
            for i in range(np.size(UD)-1):
                if (UD[i]==0 and UD[i+1]==1):
                    D2U.append(i)
                elif (UD[i]==1 and UD[i+1]==0):
                    U2D.append(i)
    
            D2U = np.asarray(D2U)
            U2D = np.asarray(U2D)  
            UDs.append(UD)
            for i in range(1,np.size(D2U)-1):
                temp = U2D[U2D>D2U[i]]   
                if (np.size(temp)>0 and temp[0] -D2U[i]>200):
                    Up = x[num,D2U[i]:temp[0]]
                    trans = []

                    #peaks = find_peaks(zscore(smooth(Up, window_len=round(200*17))), distance=160*17, height=.3)
                    peaks = find_peaks(zscore(smooth(Up, window_len=round(200))), distance=160, height=.6)
                    
                    ####
                    #This one gets a better corr with Glu, so I need to check if I'm
                    #detecting all peaks correctly
                    #peaks = find_peaks(zscore(smooth(Up, window_len=round(200*17))), distance=160*17, height=1)
                    ####
                    trans.append(np.size(peaks[0]))
                    trans.append(np.size(Up))
                    trans = np.asarray(trans)
                    features.append(np.transpose(trans))
                    cell.append(num)
    features = np.asarray(features)

    
    parameters = np.zeros((np.size(np.unique(cell)),2))
    for i in range(np.size(np.unique(cell))):                  
        for j in range(2):
                index = [k for (k,val) in enumerate(cell) if val==i]
                parameters[i,j] = np.mean(features[index,j])




    return parameters, UDs


def obtain_onsets(UDs):
    D2U  = []
    U2D = []
    for k in range(len(UDs)):
        temp_D2U = []
        temp_U2D = []
        for i in range(len(UDs[k])-1):
            if ((UDs[k][i]==0) and (UDs[k][i+1]==1)):
                temp_D2U.append(i)
            elif ((UDs[k][i]==1) and (UDs[k][i+1]==0)):
                temp_U2D.append(i)
        D2U.append(np.asarray(temp_D2U))
        U2D.append(np.asarray(temp_U2D))
                   
    return D2U, U2D
    

#%%
main_folder = r"G:\My Drive\Manuscript SWO"

meta_info = pd.read_excel(main_folder+'/LFPs information.xlsx',sheet_name='DLS')
LFPs_log = np.vstack((meta_info['Pfc'],meta_info['M1'],meta_info['BFl'],meta_info['V1']))
DLS = loadmat(main_folder+'/DLS_with_LFPs.mat')
DLS = DLS['DLS_with_LFPs']

MSN_traces = []
all_heatmaps = []
all_parameters = []
for i in range(0,len(DLS.T)):
    MSN_traces.append(DLS[0,i][-1,:])
""" 
values = []
values.append(np.mean(MSN_traces,axis=0))
values.append(np.mean(MSN_traces,axis=0)) """
for k in range(0,len(MSN_traces)):
    mua = SpectralMUA_library.computeMUA(MSN_traces[k])
    logMUA = np.log(mua)
    logMUA = SpectralMUA_library.smooth(logMUA,16) #80 ms smoothing as in Capone Cerebral Cortex 2020
    logMUA = logMUA[8:]
    parameters, Pfc_UDs_MUA = SpectralMUA_library.compute_UDs_logMUA(np.asarray(np.vstack((logMUA,logMUA))))
    #Note that the above line is extracting peaks from the MUA...
    #It can be correct, but I still need to validate it
    Pfc_UDs = []
    for i in range(0,len(Pfc_UDs_MUA)):
        Pfc_UDs[i] = np.repeat(Pfc_UDs_MUA[i],5) #Upsample UDs to 1 KHz check them in the LFP when required
        #Note that I'm upsampling to 1 KHz, but I need 20k to extract properly the 
        # MUA, so original LFPs should be downsampled to match
#%%
#######################
####OLD VERSIONS#######
#######################

for k in range(0,len(MSN_traces)):
    parameters, Pfc_UDs = compute_parameters(np.asarray(np.vstack((MSN_traces[k],MSN_traces[k]))))
    D2U_Pfc,U2D_Pfc = obtain_onsets(Pfc_UDs)
    all_parameters.append(parameters)
    psth = []
    traces = DLS[0,k]
    temp_psths = []
    for i in range(0,len(traces)):
        transitions = D2U_Pfc[0]
        psth_temp = np.zeros((len(transitions),1000))
        for k in range(2,len(transitions)-2):
            psth_temp[k,:] = traces[i,transitions[k]-500:transitions[k]+500]
        temp_psths.append(psth_temp[2:-2,:])
    all_heatmaps.append(temp_psths)

from sklearn.decomposition import PCA 

all_heatmaps_ordered = all_heatmaps[:]
for i in range(0,len(all_heatmaps_ordered)):
    neuron = all_heatmaps_ordered[i]
    for k in range(0,len(neuron)):
        test = neuron[-1]
        #test = test[2:-1,:]
        pca = PCA(n_components=5)
        pca.fit(test.T)
        #test = test[np.argsort(pca.components_[0,:]),:]
        index = np.argsort(pca.components_[0,:])

        all_heatmaps_ordered[i][k] = all_heatmaps_ordered[i][k][index,:]


#%%
plt.figure()
plt.subplot(121)
plt.imshow(all_heatmaps_ordered[0][-1],aspect = 'auto')
plt.subplot(122)
plt.imshow(all_heatmaps_ordered[0][0],aspect = 'auto')

#%% Plotting different experiments
meta_info = pd.read_excel(main_folder+'/LFPs information.xlsx',sheet_name='DLS')
cortical = np.asarray(meta_info[['Pfc','M1','BFl','V1']])

recordings = np.array([1,5,16])
subplots = [3,4,7,8]
for i in recordings:
    fields = np.where(cortical[i,:]>0)[0]

    plt.figure()    
    plt.subplot(121)
    plt.imshow(all_heatmaps_ordered[i][-1],aspect = 'auto')
    plt.title('Medium Spiny Neuron',fontsize = 25)
    plt.ylabel('Up state',fontsize = 20)
    plt.xlabel('Time (ms)',fontsize = 20)
    plt.subplot(2,4,subplots[fields[0]])
    plt.imshow(-all_heatmaps_ordered[i][0],aspect = 'auto')
 
    if len(fields)>1:
        plt.subplot(2,4,subplots[fields[1]])
        plt.imshow(-all_heatmaps_ordered[i][1],aspect = 'auto')


    plt.subplot(2,4,3),plt.title('FrA',fontsize = 25)
    plt.ylabel('Up state',fontsize = 20)

    plt.subplot(2,4,4),plt.title('M1',fontsize = 25)
    plt.subplot(2,4,7),plt.title('BF',fontsize = 25)
    plt.ylabel('Up state',fontsize = 20)
    plt.xlabel('Time (ms)',fontsize = 20)

    plt.subplot(2,4,8),plt.title('V1',fontsize = 25)
    plt.xlabel('Time (ms)',fontsize = 20)

#%% Up to Up xcorrelation
from scipy import signal
from scipy.stats import zscore
from dtaidistance import dtw

meta_info = pd.read_excel(main_folder+'/LFPs information.xlsx',sheet_name='DLS')
cortical = np.asarray(meta_info[['Pfc','M1','BFl','V1']])
Pfc_xcorr = []
M1_xcorr = []
BFl_xcorr = []
V1_xcorr = []
recordings = np.arange(0,len(meta_info))
for i in recordings:
    fields = np.where(cortical[i,:]>0)[0]
    for k in range(0,len(fields)):
        lfp = all_heatmaps_ordered[i][k]
        msn = all_heatmaps_ordered[i][-1]
        temp_corrs = np.zeros(len(msn))
        for j in range(0,len(msn)):
            temp_corrs[j] = max(signal.correlate(zscore(lfp[j,500:]),zscore(msn[j,500:]))/1000)
            #temp_corrs[j] = dtw.distance_fast(zscore(lfp[j,500:]),zscore(msn[j,500:]))

        if(fields[k]==0): Pfc_xcorr.append(np.mean(temp_corrs))
        if(fields[k]==1): M1_xcorr.append(np.mean(temp_corrs))
        if(fields[k]==2): BFl_xcorr.append(np.mean(temp_corrs))
        if(fields[k]==3): V1_xcorr.append(np.mean(temp_corrs))

plt.figure()
plt.boxplot(Pfc_xcorr,positions=[1])
plt.boxplot(M1_xcorr,positions=[2])
plt.boxplot(BFl_xcorr,positions=[3])
plt.boxplot(V1_xcorr,positions=[4])


#%% Plotting Up length relationship thourgh 2D plots


#%% Plotting different experiments
meta_info = pd.read_excel(main_folder+'/LFPs information.xlsx',sheet_name='DLS')
cortical = np.asarray(meta_info[['Pfc','M1','BFl','V1']])

recordings = np.array([1,5,16])
subplots = [1,2,3,4]
for i in recordings:
    fields = np.where(cortical[i,:]>0)[0]

    plt.figure()    
    #plt.subplot(121)
    #plt.title('Medium Spiny Neuron',fontsize = 25)
    #plt.ylabel('Up state',fontsize = 20)
    #plt.xlabel('Time (ms)',fontsize = 20)
    plt.subplot(2,2,subplots[fields[0]])
    plt.plot(all_heatmaps_ordered[i][-1][:10,:],-all_heatmaps_ordered[i][0][:10,:],'b')
    plt.plot(all_heatmaps_ordered[i][-1][-10:,:],-all_heatmaps_ordered[i][0][-10:,:],'r')
 
    if len(fields)>1:
        plt.subplot(2,2,subplots[fields[1]])
        plt.plot(all_heatmaps_ordered[i][-1][:10,:],-all_heatmaps_ordered[i][1][:10,:],'b')
        plt.plot(all_heatmaps_ordered[i][-1][-10:,:],-all_heatmaps_ordered[i][1][-10:,:],'r')

    plt.subplot(2,4,3),plt.title('FrA',fontsize = 25)
    plt.ylabel('Up state',fontsize = 20)

    plt.subplot(2,4,4),plt.title('M1',fontsize = 25)
    plt.subplot(2,4,7),plt.title('BF',fontsize = 25)
    plt.ylabel('Up state',fontsize = 20)
    plt.xlabel('Time (ms)',fontsize = 20)

    plt.subplot(2,4,8),plt.title('V1',fontsize = 25)
    plt.xlabel('Time (ms)',fontsize = 20)