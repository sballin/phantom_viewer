import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


def load_video(shot, camera):
    tree = MDSplus.Tree('spectroscopy', shot)
    return tree.getNode('gpi.%s.frames' % camera).data()


def animate_video(frames):
    fig = plt.figure()
    im = plt.imshow(frames[0], cmap=plt.get_cmap('gray'))
    
    def init():
        im.set_data(frames[0])
        return [im]
    
    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=30000, interval=0, blit=True)
    plt.show()


def flip_horizontal(frames):
    return frames[:,:,::-1]


def average_frames(frames, interval):
    frame_count = frames.shape[0]
    # Average for every [interval] frames, length = frame_count/interval
    avg_count = frame_count/interval
    avg = np.mean(frames.reshape(interval, avg_count, 64, 64), axis=0) 
    # Repeat average for every [interval] frames
    avg = np.array([[avg[i] for j in xrange(interval)] for i in 
                    xrange(avg_count)]).reshape(frame_count, 64, 64)
    return frames - avg


frames = load_video(1150528015, 'phantom2')
frames = flip_horizontal(frames)
frames = average_frames(frames, 10)
animate_video(frames)
