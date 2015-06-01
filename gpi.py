import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import numpy as np


def get_series(shot, camera, series_name):
    tree = MDSplus.Tree('spectroscopy', shot)
    if series_name == 'time': 
        try: series = tree.getNode('gpi.%s.t_hists' % camera).dim_of().data()
        except MDSplus._tdishr.TdiException: 
            sys.exit('ERROR loading time')
    else: 
        try: series = tree.getNode('gpi.%s.%s' % (camera, series_name)).data()
        except MDSplus._tdishr.TdiException: 
            sys.exit('ERROR loading %s' % series_name)
    print 'Loaded %s' % series_name
    return series 


def flip_horizontal(frames):
    return frames[:,:,::-1]


def average_frames(frames, interval):
    dim = frames.shape
    # Average for every [interval] frames, length = frame_count/interval
    avg_count = dim[0]/interval
    avg = np.mean(frames.reshape(interval, avg_count, dim[1], dim[2]), axis=0) 
    # Averages [a b c d] => [a a a... b b b... c c c... d d d...] 
    return np.array([[avg[i] for j in xrange(interval)] for i in 
                     xrange(avg_count)]).reshape(dim)


def subtract_average(frames, interval):
    # Truncate video to largest multiple of interval
    remainder = frames.shape[0] % interval 
    if remainder != 0: frames = frames[:-remainder]
    return frames - average_frames(frames, interval)


def sum_frames(frames):
    sum_frames = [np.sum(frames[i]) for i in range(frames.shape[0])]
    plt.figure()
    plt.plot(sum_frames)
    plt.show()


def animate_video(frames):
    fig = plt.figure()
    im = plt.imshow(frames[0], cmap=plt.get_cmap('gray'))
    
    def init():
        im.set_data(frames[0])
        return [im]
    
    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=0, blit=True)
    plt.show()


def slide_frames(time, frames):
    frame_count = frames.shape[0]
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[0], cmap=plt.get_cmap('gray'))
    slide_area = plt.axes([0.10, 0.1, 0.70, 0.03])
    slider = Slider(slide_area, 'Time', time[0], time[-1], valinit=time[0])
    slider.drawon = True
    slider.valfmt = '%.5f'
    curr_frame = 0

    def update(val):
        global curr_frame
        curr_frame = int((slider.val-time[0])/(time[-1]-time[0])*frame_count)
        slider.valtext.set_text('%.5f (f %d)' % (val, curr_frame))
        if curr_frame < frame_count: im.set_array(frames[curr_frame])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    playbutton_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    playbutton = Button(playbutton_area, 'Play')

    def init():
        global curr_frame
        start_frame = curr_frame
        im.set_data(frames[curr_frame])
        for i in xrange(curr_frame, curr_frame + 100): 
            slider.set_val(i/float(frame_count)*(time[-1]-time[0])+time[0])
        slider.set_val(start_frame/float(frame_count)*(time[-1]-time[0])+time[0])
        return [im]

    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=30000, interval=0, blit=True)

    playbutton.on_clicked(play)
    plt.show()


shot = 1150528015   # Cziegler: 1101209014 with apd_array 
camera = 'phantom2'
frames = get_series(shot, camera, 'frames') 

frames = flip_horizontal(frames)
frames = subtract_average(frames, 5)

time = get_series(shot, camera, 'time')
slide_frames(time, frames)

