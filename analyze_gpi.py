import sys
import MDSplus
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import signals
import eqtools
import gpi


def pixel_history(frames, x, y):
    """Return pixel x, y's values from frame array of shape (frame_count, pixels_x, pixels_y)
    """
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


# Cziegler: 1101209014 with apd_array 
shot = 1150528015   
camera = 'phantom2'

time = gpi.get_series(shot, camera, 'time')

centers = [(54, 32), (40, 50), (10, 32), (32, 10)]

for (x, y) in centers:
    pixel = np.zeros(gpi.frames.shape[0])
    region = surrounding_pixels(x, y, 5)
    for p in region:
        pixel += gpi.frames[:, p[0], p[1]] 

    before_transition = gpi.find_nearest(time, .61329)
    pixel_before = pixel[0:before_transition]
    time_before = time[0:before_transition]
    
    after_transition = gpi.find_nearest(time, .61601)
    pixel_after = pixel[0:after_transition]
    time_after = time[0:after_transition]
    
    freqs_after, PS_after = signals.power_spectrum(time_after, pixel_after)
    freqs_before, PS_before = signals.power_spectrum(time_before, pixel_before)

    plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
    plt.subplot(gs[0])
    show_region(gpi.frames, region)
    plt.subplot(gs[1])
    plt.title('Power spectrum for points around (%s, %s)' % (x, y))
    plt.plot(freqs_after, PS_after, c='r', label='after')
    plt.plot(freqs_before, PS_before, c='b', label='before')
    plt.legend()
    plt.yscale('log')
    plt.ylabel('Magnitude')
    plt.xlabel('Frequency (Hz)')
    plt.autoscale()
    plt.tight_layout(pad=1)
    plt.show()
    
