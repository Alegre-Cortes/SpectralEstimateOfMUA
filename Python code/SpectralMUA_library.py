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

#%%
import numpy as np
from scipy.signal import detrend, welch

"""Original functions from Maurizio, converted to Python using Gemini"""

def plot_median_psd_of_lfp(lfp, moving_window_size, detrending_order=2):
  """
  Plots the median power spectral density (PSD) of Local Field Potentials (LFPs).

  Args:
      lfp: A dictionary containing LFP data.
          - lfp['time']: Time vector.
          - lfp['value']: LFP values.
      moving_window_size: Window size for spectral analysis (in seconds).
      detrending_order: Order of polynomial detrending (0 to 3).

  Returns:
      F: Frequencies (Hz).
      PSDMedian: Median PSD (microvolts^2/Hz).
  """

  sampling_freq = 1 / np.diff(lfp['time'][1:2])
  window_size = round(moving_window_size * sampling_freq)
  sample_num = int(np.floor(len(lfp['value']) / window_size))

  lfps = lfp['value'][1:sample_num * window_size].reshape(window_size, sample_num)

  # Detrending
  if detrending_order > 0:
    lfps -= np.mean(lfps, axis=0, keepdims=True)  # Remove linear trend
  if detrending_order > 1:
    xx = np.linspace(-window_size // 2, window_size // 2, window_size).reshape(-1, 1)
    lfps -= xx * np.mean(np.diff(lfps, axis=0), axis=0, keepdims=True)  # Remove quadratic trend
  if detrending_order > 2:
    lfps -= 0.5 * xx**2 * np.mean(np.diff(lfps, axis=1), axis=0, keepdims=True)  # Remove cubic trend
  if detrending_order > 3:
    lfps -= 1 / 6 * xx**3 * np.mean(np.diff(lfps, axis=2), axis=0, keepdims=True)  # Remove quartic trend

  # Spectral estimation
  f, pyy = welch(lfps, fs=sampling_freq, window='hanning', nperseg=window_size)

  # Median PSD
  psd_median = np.percentile(pyy.T, 50, axis=0)

  # Plotting
  import matplotlib.pyplot as plt

  plt.figure()
  plt.hold(True)

  # Percentiles
  for prc in range(10, 50, 10):
    clr = [1 - prc / 100, 1 - prc / 100, 1]
    yp = np.percentile(pyy.T, [prc, 100 - prc], axis=0)
    plt.fill_between(f, yp[0], yp[1], color=clr, edgecolor='none')

  # Median
  plt.plot(f, psd_median, '.-b', linewidth=1)

  plt.xlim(f[0], f[-1])
  plt.yscale('log')
  plt.xscale('log')
  plt.gca().tick_params(which='both', direction='out', top=True, bottom=True)
  plt.xlabel('omega/2\pi (Hz)')
  plt.ylabel('LFP p.s.d. (uV^2/Hz)')
  plt.show()

  return f, psd_median


def compute_spectral_estimate_of_mua(t, lfp, muafreq_band, plot_power_spectrum=False, pyy_baseline=None):
  """
  Computes the spectral estimate of the Multi-Unit Activity (MUA) within a frequency band.

  Args:
      t: Time vector.
      lfp: Local Field Potential (LFP) signal.
      muafreq_band: Frequency band of interest for MUA (Hz, 2-element list).
      plot_power_spectrum: Boolean flag to plot the baseline power spectrum.
      pyy_baseline: Optional baseline power spectrum (default: None).

  Returns:
      t: Time vector downsampled to match the spectral estimate.
      MUA: Mean power spectrum of MUA within the frequency band.
  """

  dt = t[1] - t[0]
  fs = 1 / dt  # Sampling frequency
  fft_window_size = round(fs / muafreq_band[0])
  sample_num = int(np.floor(len(lfp) / fft_window_size))

  x = lfp[:fft_window_size * sample_num].reshape(fft_window_size, sample_num)

  # Detrending
  if DETRENDING_ORDER > 0:
    x -= np.mean(x, axis=0, keepdims=True)  # Remove linear trend
  if DETRENDING_ORDER > 1:
    xx = np.linspace(-fft_window_size // 2, fft_window_size // 2, fft_window_size).reshape(-1, 1)
    x -= xx * np.mean(np.diff(x, axis=0), axis=0, keepdims=True)  # Remove quadratic trend
  if DETRENING_ORDER > 2:
    x -= 0.5 * xx**2 * np.mean(np.diff(x, axis=1), axis=0, keepdims=True)  # Remove cubic trend
  if DETRENING_ORDER > 3:
    x -= 1 / 6 * xx**3 * np.mean(np.diff(x, axis=2), axis=0, keepdims=True)  # Remove quartic trend

  # Power spectrum
  f, pyy = welch(x, fs=fs, window='hanning', nperseg=fft_window_size)

  # Frequency range of interest
  ndxF = np.where((f > muafreq_band[0] - 2.0) & (f < muafreq_band[1] + 2.0))[0]

  # Baseline power spectrum
  if pyy_baseline is None:
    pyy_baseline = np.mean(pyy[ndxF, :], axis=1)
    pyy_baseline /= pyy_baseline[ndxF[0]]
    if pyy_baseline.shape[0] < pyy_baseline.shape[1]:
      pyy_baseline = pyy_baseline.T

  # Plotting
  if plot_power_spectrum:
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(f[ndxF], pyy_baseline, '.-')
    plt.xlim([f[2], f[-1]])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('omega/2\pi (Hz)')
    plt.ylabel('baseline Power spectrum (a.u.)')
    plt.show()

  # Normalized power spectrum and MUA
  norm_pyy = pyy[ndxF, :] / np.tile(pyy_baseline, (len(ndxF), 1))
  mua = np.mean(norm_pyy, axis=0)

  # Downsampled time vector
  t = np.mean(t[:fft_window_size * sample_num].reshape(fft_window_size, sample_num), axis=0)

  return t, mua
#%% 
def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    """

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y



def computeMUA(x, sampling_freq = 20000,window_size = 5, muafreq_band = [200, 1500]):
    """Extracts the MUA stimate of an individual signal following 
    the pipeline shared by Maurizio 26/02/2024"""
    window = int(window_size*sampling_freq/1000)
    x = np.reshape(x,(int(np.round(len(x)/window)),window))

    f, pyy = welch(x, fs=sampling_freq, window='hann', nperseg=window)

    # Median PSD
    psd_median = np.percentile(pyy, 50, axis=0)

    # Compute MUA
    ndxF = np.where((f > muafreq_band[0] - 2.0) & (f < muafreq_band[1] + 2.0))[0]
    norm_pyy = pyy[:,ndxF] / np.tile(psd_median[ndxF],(len(pyy),1))
    mua = np.mean(norm_pyy, axis=1)
    return mua




def compute_UDs_logMUA(dataset):
#      global dataset, parameters
    cell = []
    features = []
    UDs = []
    x = dataset
    for num in range(np.squeeze(np.shape(x[:,0]))):
        
            a = smooth(x[num,:], window_len=24)
            a = a[12:]
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
                if (np.size(temp)>0 and temp[0] - U2D[i]<50):
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
                if (np.size(temp)>0 and temp[0]-D2U[i]<4):
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
                if (np.size(temp)>0 and temp[0] -D2U[i]>40):
                    Up = x[num,D2U[i]:temp[0]]
                    trans = []

                    #peaks = find_peaks(zscore(smooth(Up, window_len=round(200*17))), distance=160*17, height=.3)
                    peaks = find_peaks(zscore(smooth(Up, window_len=round(40))), distance=32, height=.6)
                    
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