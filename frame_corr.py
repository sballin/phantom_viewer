import matplotlib.pyplot as plt
import numpy as np
import signals
import acquire
import view
import norm_xcorr


def corr_lag(t_hists, ref, time, nofft=False):
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
            xcorr = xcorr[len(xcorr)/2:len(xcorr)/2+41] 
        plt.plot(xcorr)
        plt.annotate('%s %s' % (ref, j), xy=(np.argmax(xcorr), np.max(xcorr))) 
    plt.xlabel('Lag'); plt.ylabel('Magnitude')
    plt.show()
    

def corr_frame(frames, pixel, nofft=False, other_pixels=None):
    """
    Compute the -1 to 1 normalized cross-correlation between a given pixel and
    all other pixels in the frame over the length of the video. Display a color 
    plot of the correlations for a certain lag and indicate the pixel used.
    """
    ref = frames[:, pixel[0], pixel[1]]
    #freqs, _ = signals.csd(ref, ref, return_onesided=True, detrend='linear')
    max_corrs = np.zeros((64, 64))
    if nofft: corr = np.zeros((10, 64, 64))
    else: corr = np.zeros(frames.shape)
    for i in range(64):
        for j in range(64):
            if nofft: corr[:, i, j] = [signals.cross_correlation(ref, frames[:, i, j], lag=f) for f in np.arange(0, 10)] 
            else: corr[:, i, j] = norm_xcorr.norm_xcorr(ref, frames[:, i, j]) 
            max_corrs[i, j] = np.max(corr[:, i, j])
        print i, '/ 63'
    plt.figure(); plt.imshow(max_corrs, origin='bottom'); plt.colorbar(); plt.show()
    view.slide_corr(corr, pixel, other_pixels=other_pixels)


def corr_lag_custom(frames, pixel, other_pixels):
    ref = frames[:, pixel[0], pixel[1]]
    plt.figure()
    for i, o in enumerate(other_pixels):
        signal = frames[:, o[0], o[1]]
        xcorr = [signals.cross_correlation(ref, signal, lag=f)
                 for f in np.arange(41)]
        plt.plot(xcorr)
        plt.annotate(str(i), xy=(np.argmax(xcorr), np.max(xcorr))) 
    plt.xlabel('Lag'); plt.ylabel('Magnitude')
    plt.show()


if __name__ == '__main__':
    shot = 1150611004 #1150528015 
    camera = 'phantom2'
    frames = acquire.video(shot, camera, sub=5, sobel=True)
    t_hists = acquire.gpi_series(shot, camera, 't_hists')#[:, :4]
    time = acquire.gpi_series(shot, camera, 'time')

    corr_lag_custom(frames, (4, 4), [(i, i) for i in range(5, 13)])
    corr_frame(frames[22500:24325], (8, 9), other_pixels=[(i, i) for i in range(5, 13)])
    #corr_lag(t_hists[22500:24325], 1, time, nofft=True)
 
