import sys
import MDSplus
import matplotlib.pyplot as plt
import numpy as np
import signals
import eqtools
import gpi
import scipy.signal
import bicoherence


def pixel_history(frames, x, y):
    return np.take(np.take(frames, x, axis=1), y, axis=1)


def surrounding_pixels(x, y, side):
    pixels = set()
    range = np.arange(-side/2 + 1, side/2 + 1)
    for xi in range: 
        for yi in range: pixels.add((x+xi, y+yi))
    return pixels


def show_region(frames, region):
    first_frame = np.zeros((64, 64))
    first_frame = frames[0]
    for pixel in region:
        first_frame[pixel[0]][pixel[1]] = 500
    #plt.figure()
    plt.imshow(first_frame, origin='bottom')
    plt.xlabel('pixel y coordinate')
    plt.ylabel('pixel x coordinate')
    #plt.show()


def PS_analysis(shot, camera, frames, centers):
    time = gpi.get_gpi_series(shot, camera, 'time')
    
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
        
        time_step = (time[-1]-time[0])/len(time)
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
    time = gpi.get_gpi_series(shot, camera, 'time')
    
    bicohs = []
    for (x, y) in centers:
        pixel = np.zeros(frames.shape[0])
        region = surrounding_pixels(x, y, 5)
        for p in region:
            pixel += frames[:, p[0], p[1]] 
    
        after_transition = gpi.find_nearest(time, .61601)
        pixel = pixel[after_transition:]
        time = time[after_transition:]
     
        pixel = np.array(pixel[:128*128]).reshape(32, 512) # truncate
        pixel = np.swapaxes(pixel, 0, 1)
    
        time_step = (time[-1]-time[0])/len(time)
        bicohs.append(bicoherence.bicoherence(pixel, time_step, nfft=1024, disp=False))
    
    plt.figure()
    i = 1
    for bicoh, waxis in bicohs:
        plt.subplot(2, 2, i)
        plt.title('point %d' % i)
        i += 1
        plt.xlabel('f1')
        plt.ylabel('f2')
        cont = plt.contourf(waxis, waxis, abs(bicoh), 100, cmap=plt.cm.Spectral_r)
        plt.colorbar(cont)
    
    plt.tight_layout(pad=1)
    plt.show()

 
shot = 1150528015   
camera = 'phantom2'
#frames = gpi.flip_horizontal(gpi.get_gpi_series(shot, camera, 'frames'))
centers = [(54, 32), (40, 50), (10, 32), (32, 10)]
#bicoh_analysis(shot, camera, frames, centers)
#PS_analysis(shot, camera, frames, centers)

# This is the procedure to change the number of records
qpc = np.swapaxes(scipy.io.loadmat('qpc.mat')['zmat'], 0, 1)
qpc = qpc[:32]
qpc = np.swapaxes(qpc, 0, 1)
bicoherence.bicoherence(qpc, 1, disp=True)

# This shows the power spectra for the quadratic phase coupling problem
qpc = np.swapaxes(qpc, 0, 1)
freqs, ps = scipy.signal.periodogram(qpc[1])
for i in range(qpc.shape[1]):
    ps += scipy.signal.periodogram(qpc[i])[1] 
plt.figure()
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.plot(freqs, ps)
plt.show()

