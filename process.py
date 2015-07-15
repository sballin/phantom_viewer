import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.ndimage


def time_crop((time_s, signal), time):
    """
    Crop supplied timepoints and signal to given timepoints.
    Parameters
        (time_s, signal): arrays, timepoints and signal to crop
        time: array, timepoints within which to crop
    Returns
        (time_s, signal): arrays, cropped timepoints and signal
    """
    argmin = find_nearest(time_s, time[0])
    argmax = find_nearest(time_s, time[-1])
    return time_s[argmin:argmax], signal[argmin:argmax]


def find_nearest(array, value):
    """
    Find index of first value in array closest to given value.
    Parameters
        array: NumPy array
        value: int or float
    Returns
        int: argument of array value closest to value supplied
    """
    return np.abs(array - value).argmin()


def flip_horizontal(frames):
    """
    Flip frame array horizontally.
    Parameters
        frames: frame array with dimension (frame count, y pixels, x pixels)
    """
    return frames[:, :, ::-1]


def average_frames(frames, interval):
    """
    Calculate a running average for the given frames.
    Parameters
        frames: NumPy array with dimension (frame count, y pixels, x pixels)
        interval: int, number of frames to average over
    Returns
        NumPy array: repeated averages [a a a... b b b... c c c... d d d...] 
                     for each [interval] frames supplied
    """
    dim = frames.shape
    # Average for every [interval] frames, length = frame_count/interval
    avg_count = dim[0]/interval
    avg = np.mean(frames.reshape(interval, avg_count, dim[1], dim[2]), axis=0) 
    # Averages [a b c d] => [a a a... b b b... c c c... d d d...] 
    return np.array([[avg[i] for j in xrange(interval)] for i in 
                     xrange(avg_count)]).reshape(dim)


def subtract_average(frames, interval):
    """
    Subtract running average of frames to emphasize fluctuations and
    remove background. 
    Parameters
        frames: NumPy array with dimension (frame count, y pixels, x pixels)
        interval: int, number of frames to average over
    """
    # Truncate video frame count to largest multiple of interval
    remainder = frames.shape[0] % interval 
    if remainder != 0: frames = frames[:-remainder]
    return frames - average_frames(frames, interval)


def sobel(frames):
    """
    Apply a Sobel filter to the given frames.
    """
    for i in xrange(len(frames)):
        sx = scipy.ndimage.sobel(frames[i], axis=0, mode='constant')
        sy = scipy.ndimage.sobel(frames[i], axis=1, mode='constant')
        frames[i] = np.hypot(sx, sy)
    return frames 


def kill_sobel_edges(frames):
    """
    Set edges of sobel-filtered frames to constant, average values.
    """
    frames[:, :, 0] = np.ones((frames.shape[0], frames.shape[1]))*frames.mean()
    frames[:, 0, :] = np.ones((frames.shape[0], frames.shape[1]))*frames.mean()
    return frames


def sum_frames(frames):
    """
    Sum pixel values for each frame to see total changes in intensity.
    Parameters
        frames: NumPy array with dimension (frame count, y pixels, x pixels)
    Returns
        sum_frames: NumPy array with dimension (frame count)
    """
    sum_frames = [np.sum(frames[i]) for i in range(frames.shape[0])]
    plt.figure(); plt.plot(sum_frames); plt.show()
    return sum_frames


