from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def power_spectrum(time, signal):
    signal = np.array(signal)
    PS = np.abs(np.fft.fft(signal))**2
    time_step = (time[-1]-time[0])/len(time)
    freqs = np.fft.fftfreq(signal.size, time_step)
    return freqs, PS


def plot_power_spectrum(time, signal):
    freqs, PS = power_spectrum(time, signal)
    plt.figure()
    plt.plot(freqs, PS)
    plt.show()

