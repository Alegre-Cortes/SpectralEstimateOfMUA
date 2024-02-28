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
