import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import matplotlib.gridspec as gridspec
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


def get_time_ha2(shot):
    tree = MDSplus.Tree('spectroscopy', shot)
    node = tree.getNode('\\ha_2_bright')
    return node.dim_of().data(), node.data()


def get_time_dens(shot):
    tree = MDSplus.Tree('electrons', shot)
    node = tree.getNode('tci.results.nl_04')
    return node.dim_of().data(), np.array(node.data())/.6e20


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
    # Truncate video frame count to largest multiple of interval
    remainder = frames.shape[0] % interval 
    if remainder != 0: frames = frames[:-remainder]
    return frames - average_frames(frames, interval)


def sum_frames(frames):
    sum_frames = [np.sum(frames[i]) for i in range(frames.shape[0])]
    #plt.figure()
    #plt.plot(sum_frames)
    #plt.show()
    return sum_frames


def get_frame_corners(shot, camera):
    if camera == 'phantom2':
        return [[.5095, -.429], [.5825, -.431], [.583, -.357], [.51, -.355]]
    else:
        tree = MDSplus.Tree('spectroscopy', shot)
        # Convert R, Z coordinates from centimeters to meters
        return [np.array(tree.getNode('gpi.%s.image_pos.%s_corner' % (camera, corner)).data())/100.
                for corner in ['br', 'tr', 'tl', 'bl']]


def get_extents(shot, camera):
    bl, br, tr, tl = tuple(get_frame_corners(shot, camera))
    rmin = np.mean((bl[0], tl[0]))
    rmax = np.mean((br[0], tr[0]))
    zmin = np.mean((bl[1], br[1]))
    zmax = np.mean((tl[1], tr[1]))
    return [rmin, rmax, zmin, zmax]


def find_nearest(array, value):
    return np.abs(array - value).argmin()


def time_crop((time_s, signal), time):
    argmin = find_nearest(time_s, time[0])
    argmax = find_nearest(time_s, time[-1])
    return time_s[argmin:argmax], signal[argmin:argmax]
    

def animate_video(shot, camera, time, frames, efit_tree):
    frame_count = frames.shape[0]

    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = find_nearest(efit_times, time[0])
    extents = get_extents(shot, camera)
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=extents, cmap=plt.get_cmap('gray'))
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*get_frame_corners(shot, camera)), color='r')
    plt.xlim(extents[0:2])
    plt.ylim(extents[2:4])

    time_ha2, ha2 = get_time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    time_dens, dens = time_crop(get_time_dens(shot), time)
    ax2 = plt.subplot(gs[2])#.locator_params(tight=True, nbins=6)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])

    def init():
        im.set_data(frames[0])
        efit_t_index = find_nearest(efit_times, time[0])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[0])
        vl2.set_xdata(time[0])
        return [im, l, vl1, vl2]
    
    def animate(i):
        im.set_array(frames[i])
        efit_t_index = find_nearest(efit_times, i/float(frame_count)*(time[-1]-time[0]) + time[0])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, vl1, vl2]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=0, blit=True)
    plt.tight_layout(pad=1)
    plt.show()


def slide_frames(shot, camera, time, frames, efit_tree):
    frame_count = frames.shape[0]

    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = find_nearest(efit_times, time[0])
    extents = get_extents(shot, camera)
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=extents, cmap=plt.get_cmap('gray'))
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*get_frame_corners(shot, camera)), color='r')
    plt.xlim(extents[0:2])
    plt.ylim(extents[2:4])

    time_ha2, ha2 = get_time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    time_dens, dens = time_crop(get_time_dens(shot), time)
    ax2 = plt.subplot(gs[2])#.locator_params(tight=True, nbins=6)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])

    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'Frame', 0, frame_count, valinit=0)
    slider.drawon = True
    slider.valfmt = '%d'
    curr_frame = 0

    def update(val):
        global curr_frame
        curr_frame = val 
        curr_time = time[val]
        slider.valtext.set_text('%d (t=%.5f s)' % (val, curr_time))
        if val < frame_count: 
            im.set_array(frames[val])
            efit_t_index = find_nearest(efit_times, curr_time)
            l.set_xdata(rlcfs[efit_t_index])
            l.set_ydata(zlcfs[efit_t_index])
            vl1.set_xdata(curr_time)
            vl2.set_xdata(curr_time)
        fig.canvas.draw_idle()

    slider.on_changed(update)

    playbutton_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    playbutton = Button(playbutton_area, 'Play')

    def init():
        global curr_frame
        curr_frame = int(curr_frame)
        start_frame = curr_frame
        im.set_data(frames[curr_frame])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[curr_frame])
        vl2.set_xdata(time[curr_frame])
        for i in xrange(curr_frame, curr_frame + 100, 2): 
            slider.set_val(i)
            vl1.set_xdata(time[i])
            vl2.set_xdata(time[i])
        slider.set_val(start_frame)
        return [im, l, vl1, vl2]

    def animate(i):
        im.set_array(frames[i])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, vl1, vl2]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=30000, interval=0, blit=True)

    playbutton.on_clicked(play)
    plt.show()


# Cziegler: 1101209014 with apd_array 
shot = 1150528015   
camera = 'phantom2'
frames = get_series(shot, camera, 'frames') 

efit_tree = eqtools.CModEFIT.CModEFITTree(shot)

frames = flip_horizontal(frames)
frames = subtract_average(frames, 5)

time = get_series(shot, camera, 'time')
#animate_video(shot, camera, time, frames, efit_tree)
slide_frames(shot, camera, time, frames, efit_tree)

