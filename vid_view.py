import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import matplotlib.gridspec as gridspec
import numpy as np
import signals
from scipy.misc import imsave
import os


def get_frames(shot, camera='raspi2'):
    """
    Fetch frames from the given camera for the given shot.
    Parameters
        shot: int, shot number
        camera: string e.g. 'raspi2', 'matrox2', 'matrox3'
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    try: series = tree.getNode('video_daq.%s.frames' % camera).data()
    except MDSplus._tdishr.TdiException: 
        sys.exit('ERROR loading frames')
    return series 

def animate_video(frames, cmap=plt.cm.gray, disp=True, save=False):
    """
    Animate the given frames.
    Parameters
        frames: NumPy array with dimension (frame count, x pixels, y pixels)
        cmap: colormap e.g. plt.cm.hot, plt.cm.RdBu_r
        disp: whether to display animated figure
        save: whether to save animated figure to 'out.gif'
    """
    frame_count = frames.shape[0]

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[0], origin='lower', cmap=cmap)

    def init():
        im.set_data(frames[0])
        return [im]
    
    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    plt.axis('off')
    plt.tight_layout(pad=0, rect=(0, 0, 1, 1))
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frames.shape[0], interval=0, blit=True)
    if disp: plt.show()
    if save: anim.save('out.gif', writer='imagemagick', fps=20)


def output_gif(frames):
    """ 
    Save a set of frames as out.gif by producing temporary images and putting 
    them together using ImageMagick.
    Parameters
        frames: NumPy array with dimension (frame count, x pixels, y pixels)
    """
    frames = frames[:, ::-1, :]
    for i, frame in enumerate(frames):
        imsave('outframe_%05d.png' % i, frame)
    os.system('convert outframe_*.png -layers optimize out.gif')
    os.system('rm outframe_*.png')


def main():
    # Cziegler: 1101209014 with apd_array 
    shot = int(sys.argv[1])
    frames = get_frames(shot)
    #animate_video(frames, disp=True, save=True)
    tree = MDSplus.Tree('spectroscopy', shot)
    #series = tree.getNode('video_daq.matrox3.camera_2').data()
    #output_gif(frames)


if __name__ == '__main__':
    main()

