from __future__ import division
import numpy as np
import scipy
import warnings
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


def ps_explorer(series, time_step):
    freqs, PS = scipy.signal.welch(series, nperseg=256, detrend='linear', scaling='spectrum', fs=1./time_step)
    # Initial plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(bottom=0.25)
    #p, = plt.plot(freqs, PS)
    for i in range(7, 13):
        freqs, PS = scipy.signal.welch(series, nperseg=2**i, detrend='linear', scaling='spectrum', fs=1./time_step)
        plt.plot(freqs, PS)
    plt.show()


    # Slider and button settings
    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'nperseg', 7, 12, valinit=7)
    slider.drawon = True
    slider.valfmt = '2**%d'

    def update(val):
        freqs, PS = scipy.signal.welch(series, nperseg=2**val, detrend='linear', scaling='spectrum', fs=1./time_step)
        p.set_xdata(freqs); p.set_ydata(PS)
        #slider.valtext.set_text('%d' % val)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    plt.show()


def cross_correlation(a, b, lag=0):
    a_mean = a.mean(); b_mean = b.mean()
    a_fluct = a - a_mean; b_fluct = b - b_mean
    denom = np.sqrt(np.sum(a_fluct**2)*np.sum(b_fluct**2))
    if denom == 0: return 0
    if lag == 0: return np.sum(a_fluct*b_fluct)/denom 
    elif lag < 0: return np.sum((a[-lag:] - a_mean)*(b[:lag] - b_mean))/denom
    elif lag > 0: return np.sum((a[:-lag] - a_mean)*(b[lag:] - b_mean))/denom


def PS_error(signal, nperseg=256, noverlap=None):
    """
    Get standard deviation of power spectra calculated from data in a set of
    Hanning windows with linear detrend.
    Parameters
        signal: NumPy array, signal to analyze
        nperseg: number of data points per window
        noverlap: overlap of segments, half of nperseg if not provided
    Returns
        NumPy array: standard error of the mean over Hanning windows
    """
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
    """
    Calculate power by intergrating power spectrum and using variance of signal.
    """
    pspower = np.sqrt(np.sum(PS))
    varpower = np.sqrt(np.var(scipy.signal.detrend(signal, type='linear', bp=[i*signal.size//8 for i in range(8)] + [signal.size])))
    print 'Power from PS: 1/N*sqrt(sum(PS(signal)))\t= ', pspower
    print 'Power from signal: sqrt(variance(signal))\t= ', varpower


def csd(x, y, fs=1.0, window='hanning', nperseg=256, noverlap=None, nfft=None,
        detrend='constant', return_onesided=True, scaling='density', axis=-1):
    """
    Method and helpers lifted from scipy 0.16, as yet unreleased
    """

    x = (x-np.mean(x))/(np.std(x)*len(x))
    y = (y-np.mean(y))/(np.std(y)*len(y))

    freqs, Pxy, _ = _spectral_helper(x, y, fs, window, nperseg, noverlap, nfft,
                                     detrend, return_onesided, scaling, axis,
                                     mode='psd')

    # Average over windows.
    if len(Pxy.shape) >= 2 and Pxy.size > 0:
        if Pxy.shape[-1] > 1:
            Pxy = Pxy.mean(axis=-1)
        else:
            Pxy = np.reshape(Pxy, Pxy.shape[:-1])

    return freqs, Pxy


def _spectral_helper(x, y, fs=1.0, window='hanning', nperseg=256,
                    noverlap=None, nfft=None, detrend='constant',
                    return_onesided=True, scaling='spectrum', axis=-1,
                    mode='psd'):

    if mode not in ['psd', 'complex', 'magnitude', 'angle', 'phase']:
        raise ValueError("Unknown value for mode %s, must be one of: "
                         "'default', 'psd', 'complex', "
                         "'magnitude', 'angle', 'phase'" % mode)

    # If x and y are the same object we can save ourselves some computation.
    same_data = y is x

    if not same_data and mode != 'psd':
        raise ValueError("x and y must be equal if mode is not 'psd'")

    axis = int(axis)

    # Ensure we have np.arrays, get outdtype
    x = np.asarray(x)
    if not same_data:
        y = np.asarray(y)
        outdtype = np.result_type(x,y,np.complex64)
    else:
        outdtype = np.result_type(x,np.complex64)

    if not same_data:
        # Check if we can broadcast the outer axes together
        xouter = list(x.shape)
        youter = list(y.shape)
        xouter.pop(axis)
        youter.pop(axis)
        try:
            outershape = np.broadcast(np.empty(xouter), np.empty(youter)).shape
        except ValueError:
            raise ValueError('x and y cannot be broadcast together.')

    if same_data:
        if x.size == 0:
            return np.empty(x.shape), np.empty(x.shape), np.empty(x.shape)
    else:
        if x.size == 0 or y.size == 0:
            outshape = outershape + (min([x.shape[axis], y.shape[axis]]),)
            emptyout = np.rollaxis(np.empty(outshape), -1, axis)
            return emptyout, emptyout, emptyout

    if x.ndim > 1:
        if axis != -1:
            x = np.rollaxis(x, axis, len(x.shape))
            if not same_data and y.ndim > 1:
                y = np.rollaxis(y, axis, len(y.shape))

    # Check if x and y are the same length, zero-pad if neccesary
    if not same_data:
        if x.shape[-1] != y.shape[-1]:
            if x.shape[-1] < y.shape[-1]:
                pad_shape = list(x.shape)
                pad_shape[-1] = y.shape[-1] - x.shape[-1]
                x = np.concatenate((x, np.zeros(pad_shape)), -1)
            else:
                pad_shape = list(y.shape)
                pad_shape[-1] = x.shape[-1] - y.shape[-1]
                y = np.concatenate((y, np.zeros(pad_shape)), -1)

    # X and Y are same length now, can test nperseg with either
    if x.shape[-1] < nperseg:
        warnings.warn('nperseg = {0:d}, is greater than input length = {1:d}, '
                      'using nperseg = {1:d}'.format(nperseg, x.shape[-1]))
        nperseg = x.shape[-1]

    nperseg = int(nperseg)
    if nperseg < 1:
        raise ValueError('nperseg must be a positive integer')

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg//2
    elif noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    else:
        noverlap = int(noverlap)

    # Handle detrending and window functions
    if not detrend:
        def detrend_func(d):
            return d
    elif not hasattr(detrend, '__call__'):
        def detrend_func(d):
            return scipy.signal.signaltools.detrend(d, type=detrend, axis=-1)
    elif axis != -1:
        # Wrap this function so that it receives a shape that it could
        # reasonably expect to receive.
        def detrend_func(d):
            d = np.rollaxis(d, -1, axis)
            d = detrend(d)
            return np.rollaxis(d, axis, len(d.shape))
    else:
        detrend_func = detrend

    win = scipy.signal.get_window(window, nperseg)
    
    if np.result_type(win,np.complex64) != outdtype:
        win = win.astype(outdtype)

    if mode == 'psd':
        if scaling == 'density':
            scale = 1.0 / (fs * (win*win).sum())
        elif scaling == 'spectrum':
            scale = 1.0 / win.sum()**2
        else:
            raise ValueError('Unknown scaling: %r' % scaling)
    else:
        scale = 1

    if return_onesided is True:
        if np.iscomplexobj(x):
            sides = 'twosided'
        else:
            sides = 'onesided'
            if not same_data:
                if np.iscomplexobj(y):
                    sides = 'twosided'
    else:
        sides = 'twosided'

    if sides == 'twosided':
        num_freqs = nfft
    elif sides == 'onesided':
        if nperseg % 2:
            num_freqs = (nfft + 1)//2
        else:
            num_freqs = nfft//2 + 1

    # Perform the windowed FFTs
    result = _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft)
    result = result[..., :num_freqs]
    freqs = scipy.fftpack.fftfreq(nfft, 1/fs)[:num_freqs]

    if not same_data:
        # All the same operations on the y data
        result_y = _fft_helper(y, win, detrend_func, nperseg, noverlap, nfft)
        result_y = result_y[..., :num_freqs]
        result = np.conjugate(result) * result_y
    elif mode == 'psd':
        result = np.conjugate(result) * result
    elif mode == 'magnitude':
        result = np.absolute(result)
    elif mode == 'angle' or mode == 'phase':
        result = np.angle(result)
    elif mode == 'complex':
        pass

    result *= scale
    if sides == 'onesided':
        if nfft % 2:
            result[...,1:] *= 2
        else:
            # Last point is unpaired Nyquist freq point, don't double
            result[...,1:-1] *= 2

    t = np.arange(nfft/2, x.shape[-1] - nfft/2 + 1, nfft - noverlap)/float(fs)

    if sides != 'twosided' and not nfft % 2:
        # get the last value correctly, it is negative otherwise
        freqs[-1] *= -1

    # we unwrap the phase here to handle the onesided vs. twosided case
    if mode == 'phase':
        result = np.unwrap(result, axis=-1)

    result = result.astype(outdtype)

    # All imaginary parts are zero anyways
    if same_data:
        result = result.real

    # Output is going to have new last axis for window index
    if axis != -1:
        # Specify as positive axis index
        if axis < 0:
            axis = len(result.shape)-1-axis

        # Roll frequency axis back to axis where the data came from
        result = np.rollaxis(result, -1, axis)
    else:
        # Make sure window/time index is last axis
        result = np.rollaxis(result, -1, -2)

    return freqs, result, t


def _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft):
    # Created strided array of data segments
    if nperseg == 1 and noverlap == 0:
        result = x[..., np.newaxis]
    else:
        step = nperseg - noverlap
        shape = x.shape[:-1]+((x.shape[-1]-noverlap)//step, nperseg)
        strides = x.strides[:-1]+(step*x.strides[-1], x.strides[-1])
        result = np.lib.stride_tricks.as_strided(x, shape=shape,
                                                 strides=strides)

    # Detrend each data segment individually
    result = detrend_func(result)

    # Apply window by multiplication
    result = win * result

    # Perform the fft. Acts on last axis by default. Zero-pads automatically
    result = scipy.fftpack.fft(result, n=nfft)

    return result
