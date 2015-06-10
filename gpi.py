import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import numpy as np
import signals
import eqtools


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
    print 'Loaded %s for %s on shot %s' % (series_name, camera, shot)
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
    #plt.figure()
    #plt.plot(sum_frames)
    #plt.show()
    return sum_frames


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


def get_frame_corners(shot, camera):
    if camera == 'phantom2':
        return [.5095, -.429], [.5825, -.431], [.583, -.357], [.51, -.355]
    else:
        tree = MDSplus.Tree('spectroscopy', shot)
        # Convert R, Z coordinates from centimeters to meters
        br = np.array(tree.getNode('gpi.%s.image_pos.br_corner' % camera).data())/100.
        tr = np.array(tree.getNode('gpi.%s.image_pos.tr_corner' % camera).data())/100.
        tl = np.array(tree.getNode('gpi.%s.image_pos.tl_corner' % camera).data())/100.
        bl = np.array(tree.getNode('gpi.%s.image_pos.bl_corner' % camera).data())/100.
        return bl, br, tr, tl


def slide_frames(time, frames, efit_tree):
    frame_count = frames.shape[0]

    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()

    camera = 'phantom2'
    bl, br, tr, tl = get_frame_corners(shot, camera)
    rmin = np.mean((bl[0], tl[0]))
    rmax = np.mean((br[0], tr[0]))
    zmin = np.mean((bl[1], br[1]))
    zmax = np.mean((tl[1], tr[1]))

    def find_nearest(array, value):
        return np.abs(array - value).argmin()
 
    efit_t_index = find_nearest(efit_times, time[0])
 
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[0], origin='upper', extent=[rmin, rmax, zmin, zmax], cmap=plt.get_cmap('gray'))
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*[bl, br, tr, tl]), color='r')
    plt.xlim([rmin, rmax])
    plt.ylim([zmin, zmax])

    slide_area = plt.axes([0.10, 0.1, 0.70, 0.03])
    slider = Slider(slide_area, 'Time', time[0], time[-1], valinit=time[0])
    slider.drawon = True
    slider.valfmt = '%.5f'
    curr_frame = 0

    def update(val):
        global curr_frame
        curr_frame = int((slider.val-time[0])/(time[-1]-time[0])*frame_count)
        slider.valtext.set_text('%.5f (f %d)' % (val, curr_frame))
        if curr_frame < frame_count: 
            im.set_array(frames[curr_frame])
            efit_t_index = find_nearest(efit_times, val)
            l.set_xdata(rlcfs[efit_t_index])
            l.set_ydata(zlcfs[efit_t_index])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    playbutton_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    playbutton = Button(playbutton_area, 'Play')

    def init():
        global curr_frame
        start_frame = curr_frame
        im.set_data(frames[curr_frame])
        l.set_xdata(rlcfs[argt])
        l.set_ydata(zlcfs[argt])
        for i in xrange(curr_frame, curr_frame + 100): 
            slider.set_val(i/float(frame_count)*(time[-1]-time[0])+time[0])
        slider.set_val(start_frame/float(frame_count)*(time[-1]-time[0])+time[0])
        return [im, l]

    def animate(i):
        im.set_array(frames[i])
        l.set_xdata(rlcfs[argt])
        l.set_ydata(zlcfs[argt])
        return [im, l]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=30000, interval=0, blit=True)

    playbutton.on_clicked(play)
    plt.show()


def rotate(points, origin, angle):
    origin = np.array([origin[0], origin[1]])
    angle = angle*np.pi/180.
    return np.dot(points - origin, np.array([[np.cos(angle), np.sin(angle)],
                                             [-np.sin(angle), np.cos(angle)]])) + origin


# Cziegler: 1101209014 with apd_array 
shot = 1150528015   
camera = 'phantom2'
frames = get_series(shot, camera, 'frames') 

efit_tree = eqtools.CModEFIT.CModEFITTree(shot)

frames = flip_horizontal(frames)
frames = subtract_average(frames, 5)

time = get_series(shot, camera, 'time')
slide_frames(time, frames, efit_tree)

