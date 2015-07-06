import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import signals
import eqtools
import gpi
import scipy
import bicoherence
import norm_xcorr


def pixel_history(frames, x, y):
    """
    Return time history for a specified pixel.
    Parameters
        frames: NumPy array of dimension (frame count, y pixels, x pixels)
    """
    return np.take(np.take(frames, x, axis=1), y, axis=1)


def surrounding_pixels(x, y, side):
    """
    Return array of (x, y) pairs surrounding pixel at (x, y).
    Parameters
        side: length of box around pixel
    """
    pixels = set()
    range = np.arange(-side/2 + 1, side/2 + 1)
    for xi in range: 
        for yi in range: pixels.add((x+xi, y+yi))
    return pixels


def show_region(frames, region):
    """
    Display set of pixels used on top of the first frame.
    Parameters
        frames: NumPy array of dimension (frame count, y pixels, x pixels)
        region: array of (x, y) tuples specifying used pixels
    """
    first_frame = np.zeros((frames.shape[1], frames.shape[2]))
    first_frame = frames[0].astype(float)
    for pixel in region:
        first_frame[pixel[0], pixel[1]] = np.nan
    #plt.figure()
    cmap = plt.cm.gray
    cmap.set_bad((1, 0, 0, 1))
    plt.imshow(first_frame, origin='bottom', cmap=cmap, interpolation='nearest')
    plt.xlabel('pixel y coordinate')
    plt.ylabel('pixel x coordinate')
    #plt.show()


def PS_analysis(shot, camera, frames, centers):
    """
    Display signal and power series for centers specified.
    """
    time = gpi.get_gpi_series(shot, camera, 'time')
    time_step = (time[-1]-time[0])/len(time)
    print 1./time_step
    
    for (x, y) in centers:
        pixel = np.zeros(frames.shape[0])
        region = surrounding_pixels(x, y, 5)
        for p in region:
            pixel += frames[:, p[0], p[1]] 
    
        winlen = 1024
        freqs, PS = scipy.signal.welch(pixel, fs=1./time_step, nperseg=winlen, detrend='linear', scaling='spectrum')

        print 'Point', x, y
        print 'FFT window length: %d' % winlen
        signals.power_analysis(pixel, PS)

        plt.figure()
        plt.subplot2grid((2, 2), (0,0))
        show_region(frames, region)
        plt.subplot2grid((2, 2), (0, 1))
        plt.plot(time, pixel)
        plt.xlabel('time')
        plt.ylabel('signal')
        plt.autoscale(tight=True)
        plt.subplot2grid((2, 2), (1, 0), colspan=2)
        plt.title('Power spectrum for points around (%s, %s)' % (x, y))
        plt.semilogy(freqs, PS, 'b-')
        #plt.xscale('log')
        error = signals.PS_error(pixel, nperseg=winlen)
        plt.fill_between(freqs, PS-error, PS+error, color='b', alpha=.5)
        plt.ylabel('Magnitude')
        plt.xlabel('Frequency (Hz)')
        plt.autoscale(tight=True)
        plt.tight_layout(pad=1)
        
    plt.show()
    

def split_PS_analysis(shot, camera, frames, centers):
    """
    Compare power spectra before and after a certain time for pixels in the
    given centers.
    """
    time = gpi.get_gpi_series(shot, camera, 'time')
    time_step = (time[-1]-time[0])/len(time)
    
    for (x, y) in centers:
        pixel = np.zeros(frames.shape[0])
        region = surrounding_pixels(x, y, 5)
        for p in region:
            pixel += frames[:, p[0], p[1]] 
    
        before_transition = gpi.find_nearest(time, .61329)
        pixel_before = pixel[:before_transition]
        time_before = time[:before_transition]
        
        after_transition = gpi.find_nearest(time, .61601)
        pixel_after = pixel[after_transition:after_transition+time_before.size]
        time_after = time[after_transition:after_transition+time_before.size]
        
        segs = 24
        freqs_before, PS_before = scipy.signal.welch(pixel_before, fs=1./time_step, nperseg=bicoherence.nextpow2(pixel_before.size/segs), detrend='linear', scaling='spectrum')
        freqs_after, PS_after = scipy.signal.welch(pixel_after, fs=1./time_step, nperseg=bicoherence.nextpow2(pixel_after.size/segs), detrend='linear', scaling='spectrum')

        print 'Point', x, y
        signals.power_analysis(pixel_before, PS_before)
        signals.power_analysis(pixel_after, PS_after)
    
        plt.figure()
        plt.subplot2grid((2, 2), (0,0))
        show_region(frames, region)
        plt.subplot2grid((2, 2), (0, 1))
        plt.plot(time, pixel)
        plt.xlabel('time')
        plt.ylabel('signal')
        plt.subplot2grid((2, 2), (1, 0), colspan=2)
        plt.title('Power spectrum for points around (%s, %s)' % (x, y))
        plt.semilogy(freqs_after, PS_after, 'r-', label='after')
        plt.semilogy(freqs_before, PS_before, 'b-', label='before')
        error_before = signals.PS_error(pixel_before, nperseg=bicoherence.nextpow2(pixel_before.size/segs))
        error_after = signals.PS_error(pixel_after, nperseg=bicoherence.nextpow2(pixel_after.size/segs))
        plt.fill_between(freqs_before, PS_before-error_before, PS_before+error_before, color='b', alpha=.5)
        plt.fill_between(freqs_after, PS_after-error_after, PS_after+error_after, color='r', alpha=.5)
        plt.legend()
        plt.ylabel('Magnitude')
        plt.xlabel('Frequency (Hz)')
        plt.autoscale()
        plt.tight_layout(pad=1)
        
    plt.show()
    

def bicoh_analysis(shot, camera, frames, centers):
    """
    Perform a bicoherence analysis of pixels around the given centers.
    """
    time = gpi.get_gpi_series(shot, camera, 'time')
    time_step = (time[-1]-time[0])/len(time)
    
    bicohs = []
    for (x, y) in centers:
        pixel = np.zeros(frames.shape[0])
        region = surrounding_pixels(x, y, 5)
        for p in region:
            pixel += frames[:, p[0], p[1]] 
        pixel = pixel/float(len(region))
     
        pixel = np.array(pixel[:128*64]).reshape(128, 64) # truncate
        pixel = np.swapaxes(pixel, 0, 1)
    
        bicohs.append(bicoherence.bicoherence(pixel, time_step, nfft=128, disp=False))
    
    plt.figure()
    i = 1
    for bicoh, waxis in bicohs:
        plt.subplot(2, 2, i)
        plt.title('point %d' % i)
        i += 1
        plt.xlabel('f1 (kHz)')
        plt.ylabel('f2 (kHz)')
        waxis = waxis/1000.
        cont = plt.contourf(waxis, waxis, abs(bicoh), 100, cmap=plt.cm.Spectral_r, vmin=0, vmax=1)
        #plt.xlim(0, waxis[-1])
        #plt.ylim(0, waxis[-1])
        plt.colorbar(cont)
    
    plt.tight_layout(pad=1)
    plt.show()


def edge_filter(frames):
    """
    Apply a Sobel filter to the given frames.
    """
    for i in xrange(len(frames)):
        sx = scipy.ndimage.sobel(frames[i], axis=0, mode='constant')
        sy = scipy.ndimage.sobel(frames[i], axis=1, mode='constant')
        frames[i] = np.hypot(sx, sy)
    return frames 


def corr_lag(t_hists, time):
    """
    For the given pixel time histories, plot correlations as a function of
    frame lag.
    """
    t_hists = t_hists.swapaxes(0, 1)
    frame_count = t_hists.shape[1]
    plt.figure()
    for i in range(len(t_hists)):
        for j in range(len(t_hists)):
            xcorr = norm_xcorr.norm_xcorr(t_hists[i], t_hists[j])
            if i == j: plt.plot(np.arange(-frame_count/2, frame_count/2), xcorr, '--')
            else: plt.plot(np.arange(-frame_count/2, frame_count/2), xcorr)
    plt.xlabel('Lag')
    plt.ylabel('Magnitude')
    plt.show()
    

def corr_frame(frames, pixel):
    """
    Compute the -1 to 1 normalized cross-correlation between a given pixel and all 
    other pixels in the frame over the length of the video. Display a color plot of
    the correlations for a certain lag and indicate the pixel used.
    """
    ref = frames[:, pixel[0], pixel[1]]
    #freqs, _ = signals.csd(ref, ref, return_onesided=True, detrend='linear')
    corr = np.zeros(frames.shape)
    for i in range(64):
        for j in range(64):
            corr[:, i, j] = norm_xcorr.norm_xcorr(ref, frames[:, i, j]) 
        print i
    gpi.slide_corr(corr, pixel)


if __name__ == '__main__':
    shot = 1150611004 #1150528015  #1150611004 
    camera = 'phantom2'
    #frames = gpi.flip_horizontal(gpi.get_gpi_series(shot, camera, 'frames'))
    centers = [(54, 32), (40, 50), (10, 32), (32, 10)]
    
    #eframes = edge_filter(subs)
    #subs = gpi.subtract_average(frames, 5)
    #corr_frame(subs[22500:24325], (30, 30))
    
    #PS_analysis(shot, camera, frames, centers)
    t_hists = gpi.get_gpi_series(shot, camera, 't_hists')
    time = gpi.get_gpi_series(shot, camera, 'time')
    time_step = (time[-1]-time[0])/len(time)
    #corr_lag(t_hists[22500:24325], time)
    
    #signal = np.array(signal[:256*128]).reshape(128, 256) # truncate
    #signal = np.swapaxes(signal, 0, 1)
    #bicoherence.bicoherence(signal, 1., nfft=256, disp=False)
    
    #bicoh_analysis(shot, camera, frames, centers)
    
    time = gpi.get_gpi_series(shot, camera, 'time')
    before_transition = gpi.find_nearest(time, .61329)
    frames_before = frames[:before_transition]
    #bicoh_analysis(shot, camera, frames_before, centers)
    after_transition = gpi.find_nearest(time, .61601)
    frames_after = frames[after_transition:after_transition+frames_before.size]
    #bicoh_analysis(shot, camera, frames_after, centers)
     
