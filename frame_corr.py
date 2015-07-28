import matplotlib.pyplot as plt
import numpy as np
import signals
import acquire
import view
import norm_xcorr
from tqdm import *
import scipy.signal


def corr_frame(frames, pixel, nofft=False, other_pixels=None, show_maxes=False):
    """
    Compute the -1 to 1 normalized cross-correlation between a given pixel and
    all other pixels in the frame over the length of the video. Display a color 
    plot of the correlations for a certain lag and indicate the pixel used.
    """
    ref = frames[:, pixel[0], pixel[1]]
    #freqs, _ = signals.csd(ref, ref, return_onesided=True, detrend='linear')
    max_corrs = np.zeros((64, 64))
    max_lags = np.zeros((64, 64))
    #if nofft: corr = np.zeros((100, 64, 64))
    corr = np.zeros(frames.shape)
    for i in tqdm(range(64)):
        for j in range(64):
            if nofft: corr[:, i, j] = [signals.cross_correlation(ref, frames[:, i, j], lag=f) for f in np.arange(-912, 913)] 
            else: corr[:, i, j] = norm_xcorr.norm_xcorr(ref, frames[:, i, j]) 
            max_corrs[i, j] = np.max(corr[:, i, j])
            max_lags[i, j] = np.argmax(corr[:, i, j])-912

    if show_maxes:
        plt.figure(); plt.title('Maximum correlations')
        plt.imshow(max_corrs, origin='bottom'); plt.colorbar(label='Magnitude')
        plt.show()
    
        plt.figure(); plt.title('Maximum correlations')
        plt.imshow(max_lags, origin='bottom', vmin=-75, vmax=75); plt.colorbar(label='Lag (frames)')
        plt.show()

    view.slide_corr(corr, pixel, other_pixels=other_pixels)


def corr_plot(t_hists, ref, time, nofft=False):
    """
    For the given pixel time histories, plot correlations as a function of
    frame lag.
    """
    t_hists = t_hists.swapaxes(0, 1)
    frame_count = t_hists.shape[1]
    plt.figure()
    for j in range(len(t_hists)):
        if nofft: xcorr = [signals.cross_correlation(t_hists[ref], t_hists[j], 
                           lag=f) for f in np.arange(41)] 
        else: 
            xcorr = norm_xcorr.norm_xcorr(t_hists[ref], t_hists[j]) 
        plt.plot(xcorr)
        plt.annotate('%s %s' % (ref, j), xy=(np.argmax(xcorr), np.max(xcorr))) 
    plt.xlabel('Lag'); plt.ylabel('Magnitude')
    plt.show()
    

def corr_plot_custom(frames, pixel, other_pixels):
    max_lag = 500
    ref = frames[:, pixel[0], pixel[1]]
    plt.figure()
    for i in trange(len(other_pixels)):
        signal = frames[:, other_pixels[i][0], other_pixels[i][1]]
        xcorr = [signals.cross_correlation(ref, signal, lag=f)
                 for f in np.arange(-max_lag, max_lag)]
        plt.plot(np.arange(-max_lag, max_lag), xcorr)
        plt.axvline()
        plt.annotate(str(i), xy=(np.argmax(xcorr)-max_lag, np.max(xcorr))) 
    plt.xlabel('Lag'); plt.ylabel('Correlation')
    plt.show()


if __name__ == '__main__':
    shot = 1150722015 #1150611004 #1150528015 
    camera = 'phantom2'

    frames = acquire.video(shot, camera, sub=0, sobel=False)
    t_hists = acquire.gpi_series(shot, camera, 't_hists')#[:, :4]
    time = acquire.gpi_series(shot, camera, 'time')

    # Pixels along divider
    #custom_pixels = [(15, 35), (25, 34), (35, 30), (44, 24), (51, 20), (58, 13)] 
    # Pixels along blob field line
    #custom_pixels = [(15, 35), (23, 31), (30, 26), (33, 21), (36, 14), (37, 8), (35, 2)]
    # Pixels descending vertically in middle
    #custom_pixels = [(32, i) for i in reversed(xrange(4, 29, 3))]

    #corr_plot_custom(frames, (35, 5), custom_pixels)

    #corr_plot(t_hists[22500:24325], 1, time, nofft=True)

    time_step = (time[-1]-time[0])/len(time)
    plt.figure(); plt.plot(time, frames[:, 14, 10]); plt.show()
    signals.ps_explorer(frames[:, 14, 10], time_step)

    #corr_frame(frames[22500:24325], (35, 5), other_pixels=custom_pixels)
 
