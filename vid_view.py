import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import matplotlib.gridspec as gridspec
import numpy as np
import signals


def get_frames(shot):
    tree = MDSplus.Tree('spectroscopy', shot)
    try: series = tree.getNode('video_daq.raspi2.frames').data()
    except MDSplus._tdishr.TdiException: 
        sys.exit('ERROR loading frames')
    return series 

def animate_video(frames):
    frame_count = frames.shape[0]

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[0], origin='upper', cmap=plt.get_cmap('hot'))

    def init():
        im.set_data(frames[0])
        return [im]
    
    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=0, blit=True)

    plt.tight_layout(pad=1)
    plt.show()


# Cziegler: 1101209014 with apd_array 
#shot = int(sys.argv[1])
#frames = get_frames(shot)
#animate_video(frames)
tree = MDSplus.Tree('spectroscopy', 1150528015)
series = tree.getNode('video_daq.matrox3.camera_2').data()
animate_video(series)
