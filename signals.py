from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy


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

