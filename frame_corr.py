import matplotlib.pyplot as plt
import numpy as np
import signal
import acquire
import view
import norm_xcorr


def corr_lag(t_hists, time, nofft=False):
    """
    For the given pixel time histories, plot correlations as a function of
    frame lag.
    """
    t_hists = t_hists.swapaxes(0, 1)
    frame_count = t_hists.shape[1]
    plt.figure()
    for i in range(len(t_hists)):
        for j in range(len(t_hists)):
            if nofft: xcorr = [signal.cross_correlation(t_hists[i], t_hists[j], 
                               lag=f) for f in np.arange(0, 100)] 
            else: xcorr = norm_xcorr.norm_xcorr(t_hists[i], t_hists[j]) 
            if i < j: 
                new = xcorr[len(xcorr)/2:len(xcorr)/2+41] 
                plt.plot(np.arange(0, 41), new)
                plt.annotate('%s %s' % (i, j), xy=(np.argmax(new), np.max(new))) 
                #plt.plot(np.arange(-frame_count/2, frame_count/2), xcorr)
    plt.xlabel('Lag')
    plt.ylabel('Magnitude')
    plt.show()
    

def corr_frame(frames, pixel, nofft=False):
    """
    Compute the -1 to 1 normalized cross-correlation between a given pixel and
    all other pixels in the frame over the length of the video. Display a color 
    plot of the correlations for a certain lag and indicate the pixel used.
    """
    ref = frames[:, pixel[0], pixel[1]]
    #freqs, _ = signal.csd(ref, ref, return_onesided=True, detrend='linear')
    if nofft: corr = np.zeros((10, 64, 64))
    else: corr = np.zeros(frames.shape)
    for i in range(64):
        for j in range(64):
            if nofft: corr[:, i, j] = [signal.cross_correlation(ref, frames[:, i, j], lag=f) for f in np.arange(0, 10)] 
            else: corr[:, i, j] = norm_xcorr.norm_xcorr(ref, frames[:, i, j]) 
        print i
    view.slide_corr(corr, pixel)

if __name__ == '__main__':
    shot = 1150611004 #1150528015 
    camera = 'phantom2'
    frames = acquire.video(shot, camera, sub=5, sobel=True)
    t_hists = acquire.gpi_series(shot, camera, 't_hists')[:, :4]
    time = acquire.gpi_series(shot, camera, 'time')

    corr_frame(frames[22500:24325], (8, 9))
    corr_lag(t_hists[22500:24325], time)
 
