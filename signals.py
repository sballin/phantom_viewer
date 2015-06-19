from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy


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
        s -= np.mean(s)
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


def PS_error(signal, nperseg=256, noverlap=None):
    x = np.arange(signal.size)
    if not noverlap: noverlap = nperseg//2
    step = nperseg - noverlap
    shape = x.shape[:-1]+((x.shape[-1]-noverlap)//step, nperseg)
    strides = x.strides[:-1]+(step*x.strides[-1], x.strides[-1])
    indices = np.lib.stride_tricks.as_strided(x, shape=shape,
                                              strides=strides) 
    spectra = [scipy.signal.periodogram(np.hanning(indices.shape[1])
                                        *signal[indices[i]], 
                                        detrend='linear', scaling='spectrum')[1]
               for i in range(indices.shape[0])]
    return scipy.stats.sem(spectra)


def power_analysis(signal, PS):
    pspower = np.sqrt(np.sum(PS))
    varpower = np.sqrt(np.var(scipy.signal.detrend(signal, type='linear', bp=[i*signal.size//8 for i in range(8)] + [signal.size])))
    print 'Power from PS: 1/N*sqrt(sum(PS(signal)))\t= ', pspower
    print 'Power from signal: sqrt(variance(signal))\t= ', varpower

