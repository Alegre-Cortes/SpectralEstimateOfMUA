
def compute_spectral_estimate_of_mua(t, lfp, muafreqband, plot_power_spectrum=False, pyy_baseline=None):
    """
    Estimates the power spectrum of the MUA (Multi-Unit Activity) band in the LFP (Local Field Potential) signal.

    Args:
        t (np.ndarray): Time vector of the LFP signal.
        lfp (np.ndarray): LFP signal.
        muafreqband (list[float]): Frequency band of interest for MUA (e.g., [200, 1500]).
        plot_power_spectrum (bool, optional): If True, plots the baseline power spectrum. Defaults to False.
        pyy_baseline (np.ndarray, optional): Baseline power spectrum. If None, it will be calculated. Defaults to None.

    Returns:
        tuple:
            - t (np.ndarray): Time vector for the averaged MUA spectrum.
            - mua (np.ndarray): Averaged power spectrum of the MUA band.
    """

    dt = t[1] - t[0]
    fs = 1 / dt  # Sampling frequency

    # Window size based on the lower frequency bound
    fft_window_size = round(fs / muafreqband[0])
    num_windows = int(np.floor(len(lfp) / fft_window_size))

    # Reshape LFP into windows
    x = lfp[:fft_window_size * num_windows].reshape(fft_window_size, num_windows)

    # Detrending based on DETRENDING_ORDER (modify if needed)
    if 2 <= DETRENDING_ORDER <= 4:
        # Linear detrending (order 1) is assumed as the default
        x_trend = np.repeat(np.linspace(-fft_window_size // 2, fft_window_size // 2, fft_window_size), num_windows, axis=1)
        x -= x_trend * np.repeat(np.mean(np.diff(x, axis=0), axis=1), fft_window_size, axis=1)
    elif DETRENDING_ORDER > 4:
        raise ValueError("DETRENDING_ORDER must be between 0 and 4")

    # Power spectrum calculation
    y = np.fft.fft(x, axis=0)
    pyy = np.abs(y) ** 2 / fft_window_size
    f = fs * np.linspace(0, 0.5, int(fft_window_size / 2) + 1)

    # Frequency indices within the MUA band (with a slight buffer)
    f_idx = np.where((f > muafreqband[0] - 2.0) & (f < muafreqband[1] + 2.0))[0]

    # Baseline calculation
    if pyy_baseline is None:
        pyy_baseline = np.mean(pyy[f_idx, :], axis=1)
        pyy_baseline /= pyy_baseline[f_idx[0]]
    if pyy_baseline.shape[0] < pyy_baseline.shape[1]:
        pyy_baseline = pyy_baseline.T

    # Plot power spectrum if requested
    if plot_power_spectrum:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(f[f_idx], pyy_baseline, ".-")
        plt.xlim([f[2], f[-1]])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Baseline power spectrum (a.u.)")
        plt.show()

    # Normalized power spectrum and averaging
    norm_pyy = pyy[f_idx, :] / np.repeat(pyy_baseline, pyy.shape[1], axis=0)
    mua = np.mean(norm_pyy, axis=1)

    # Average time point (consider adjusting based on your needs)
    t = np.mean(t[1:fft_window_size * num_windows].reshape(fft_window_size, num_windows), axis=0)

    return t, mua