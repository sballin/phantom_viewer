from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def power_spectrum(time, signal):
    signal = np.array(signal)
    PS = np.abs(np.fft.fft(signal))**2
    time_step = (time[-1]-time[0])/len(time)
    freqs = np.fft.fftfreq(signal.size, time_step)
    idx = np.argsort(freqs)
    freqs = freqs[idx]
    PS = PS[idx]
    return freqs[len(freqs)/2:], PS[len(PS)/2:]


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

