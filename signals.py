from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt


def window(signal, winlen, overlap):
    """Apply Hanning window to a signal divided into chunks of same length.

    Args: 
        signal: signal to which windowing is applied.
        winlen: length of windows.
        overlap: overlap of segments. 0.5 = 50% overlap.

    Returns:
        Sequential overlapping signal chunks of same length.
    """
    siglen = signal.shape[0]
    han = np.hanning(winlen)
    chunks = []
    for i in np.arange(0, siglen, (1 - overlap)*winlen):
        sig = signal[i:i + winlen]
        if han.shape == sig.shape: chunks.append(han*sig)
    return chunks


def power_spectrum(time, signal):
    signal = np.array(signal)
    time = np.array(time)
    if time.shape != signal.shape or signal.shape[0] == 0: 
        sys.exit("ERROR: FFT dimension problem.")
    PS = np.abs(np.fft.fft(signal))**2
    time_step = (time[-1]-time[0])/len(time)
    freqs = np.fft.fftfreq(signal.size, time_step)
    idx = np.argsort(freqs)
    freqs = freqs[idx]
    PS = PS[idx]
    return freqs[len(freqs)/2:], PS[len(PS)/2:]


def windowed_power_spectrum(time, signal, winlen, overlap):
    time_step = (time[-1]-time[0])/len(time)
    PS_sum = 0
    for s in window(signal, winlen, overlap):
        PS_sum += np.fft.fft(s)
        freqs = np.fft.fftfreq(s.size, time_step)
    PS_sum = np.abs(PS_sum)**2
    idx = np.argsort(freqs[:winlen])
    freqs = freqs[idx]
    PS_sum = PS_sum[idx]
    return freqs[len(freqs)/2:], PS_sum[len(PS_sum)/2:]


def plot_power_spectrum(time, signal, title):
    freqs, PS = power_spectrum(time, signal)
    plt.figure()
    plt.title(title)
    plt.plot(freqs, PS)
    plt.yscale('log')
    plt.ylabel('Magnitude')
    plt.xlabel('Frequency (Hz)')
    plt.autoscale()
    plt.show()

