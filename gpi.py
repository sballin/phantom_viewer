import sys
import MDSplus
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import matplotlib.gridspec as gridspec
import numpy as np
import signals
import eqtools
import scipy


def get_gpi_series(shot, camera, series_name):
    """
    Get specified GPI-related series.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplane GPI),
                'phantom2' (X-point GPI)
        series_name: string e.g. 'time' (timepoints), 'frames' (array of frames 
                     with dimension (frame count, y pixels, x pixels)
    Returns
        series: chosen with series_name
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    if series_name == 'time': 
        try: series = tree.getNode('gpi.%s.t_hists' % camera).dim_of().data()
        except MDSplus._tdishr.TdiException: 
            print 'Time not loaded from t_hists'
            start = tree.getNode('gpi.%s.settings.trig_time' % camera).data()
            frame_rate = tree.getNode('gpi.%s.settings.frame_rate' % camera).data()
            num_frames = tree.getNode('gpi.%s.settings.num_frames' % camera).data()
            return np.arange(start, start + num_frames/float(frame_rate), 1./frame_rate)
    else: 
        try: series = tree.getNode('gpi.%s.%s' % (camera, series_name)).data()
        except MDSplus._tdishr.TdiException: 
            sys.exit('ERROR loading %s' % series_name)
    print 'Loaded %s for %s on shot %s' % (series_name, camera, shot)
    return series 


def get_time_ha2(shot):
    """
    Get timepoints and H_alpha signal for specified shot.
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    node = tree.getNode('\\ha_2_bright')
    return node.dim_of().data(), node.data()


def get_time_dens(shot):
    """
    Get timepoints and line average density for specified shot.
    """
    tree = MDSplus.Tree('electrons', shot)
    node = tree.getNode('tci.results.nl_04')
    return node.dim_of().data(), np.array(node.data())/.6e20


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


def get_frame_corners(shot, camera):
    """
    Return (R, Z) FOV corner locations for given camera in meters.
    """
    if camera == 'phantom2':
        return [[.5095, -.429], [.5825, -.431], [.583, -.357], [.51, -.355]]
    else:
        tree = MDSplus.Tree('spectroscopy', shot)
        # Convert R, Z coordinates from centimeters to meters
        return [np.array(tree.getNode('gpi.%s.image_pos.%s_corner' 
                                      % (camera, corner)).data())/100.
                for corner in ['br', 'tr', 'tl', 'bl']]


def get_extents(shot, camera):
    """
    Calculate R, Z border locations and return in matplotlib 'extent' order.
    """
    bl, br, tr, tl = tuple(get_frame_corners(shot, camera))
    rmin = np.mean((bl[0], tl[0]))
    rmax = np.mean((br[0], tr[0]))
    zmin = np.mean((bl[1], br[1]))
    zmax = np.mean((tl[1], tr[1]))
    return [rmin, rmax, zmin, zmax]


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


def get_psi_contours(shot):
    """
    Get contours of magnetic flux interpolated over a grid.
    Returns
        psinew: psi values over new grid,
        extent: (R, Z) coordinates of grid corners
    """
    tree = MDSplus.Tree('analysis', shot)
    psirz = tree.getNode('efit.results.g_eqdsk.psirz').data()
    rgrid = tree.getNode('efit.results.g_eqdsk.rgrid').data()
    zgrid = tree.getNode('efit.results.g_eqdsk.zgrid').data()
    rnew = np.arange(np.min(rgrid), np.max(rgrid), .01)
    znew = np.arange(np.min(zgrid), np.max(zgrid), .01)
    psinew = np.zeros((len(psirz), len(znew), len(rnew)))
    for i in range(len(psirz)):
        f = scipy.interpolate.interp2d(rgrid, zgrid, psirz[i], kind='cubic')
        psinew[i, :, :] = f(rnew, znew) 
    extent = [np.min(rgrid), np.max(rgrid), np.min(zgrid), np.max(zgrid)]
    return psinew, extent


def animate_video(shot, camera, time, frames):
    """
    Animate frames while displaying last closed flux surface and relevant plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        time: NumPy array of time points
        frames: NumPy array of frames with dimension (frame count, x pixels, 
                y pixels)
    """
    frame_count = frames.shape[0]

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = find_nearest(efit_times, time[0])
    gpi_extent = get_extents(shot, camera)
    psirz, psiext = get_psi_contours(shot)

    # GPI, LCFS, and flux contour initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    ax0 = plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=gpi_extent, cmap=plt.cm.gray)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    c, = plt.plot(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz), np.max(psirz), .001), extent=psiext))
    print im, l, c
    plt.scatter(*zip(*get_frame_corners(shot, camera)), color='r')
    plt.xlim(gpi_extent[0:2])
    plt.ylim(gpi_extent[2:4])

    # H_alpha timeseries plot
    time_ha2, ha2 = get_time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    # Line average density plot
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
        #c.set_array(psirz[efit_t_index])
        c.set_data(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz[0]), np.max(psirz[0]), .001), extent=psiext))
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[0])
        vl2.set_xdata(time[0])
        return [im, l, c, vl1, vl2]
    
    def animate(i):
        im.set_array(frames[i])
        efit_t_index = find_nearest(efit_times, i/float(frame_count)*(time[-1]-time[0]) + time[0])
        #c.set_array(psirz[efit_t_index])
        c.set_data(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz[0]), np.max(psirz[0]), .001), extent=psiext))
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, c, vl1, vl2]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=0, blit=True)
    plt.tight_layout(pad=1)
    plt.show()


def slide_corr(frames, pixel):
    """
    Display correlation values for a pixel with a slider for lag time.
    Parameters
        frames: NumPy array of correlations with dimension (no. lags, x pixels,
                y pixels)
        pixel: (x, y) to mark with a circle
    """
    frame_count = frames.shape[0]

    # Initial plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[0], origin='bottom', cmap=plt.cm.RdBu_r, vmin=-1, vmax=1)
    circ = plt.Circle(pixel[::-1], radius=1, edgecolor='r', fill=False)
    ax.add_patch(circ)
    plt.colorbar()
    plt.title('Correlation for pixel (%d, %d)' % pixel)
    plt.xlabel('Pixel y coordinate')
    plt.ylabel('Pixel x coordinate')

    # Slider and button settings
    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'Frame', 0, frame_count, valinit=0)
    slider.drawon = True
    slider.valfmt = '%d'
    playbutton_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    playbutton = Button(playbutton_area, 'Play')

    curr_frame = 0

    def update(val):
        global curr_frame
        curr_frame = val 
        slider.valtext.set_text('%d' % val)
        if val < frame_count: 
            im.set_array(frames[val])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    def init():
        global curr_frame
        curr_frame = int(curr_frame)
        start_frame = curr_frame
        im.set_data(frames[curr_frame])
        for i in xrange(curr_frame, curr_frame + 100, 2): 
            slider.set_val(i)
        return [im]

    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=frame_count, interval=0, blit=False)

    playbutton.on_clicked(play)
    plt.show()


def slide_frames(shot, camera, time, frames):
    """
    Slide through GPI frames while displaying last closed flux surface and 
    relevant timeseries plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        time: NumPy array of time points
        frames: NumPy array of frames with dimension (frame count, x pixels, 
                y pixels)
    """
    frame_count = frames.shape[0]

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = find_nearest(efit_times, time[0])
    extents = get_extents(shot, camera)
    
    # GPI and LCFS initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=extents, cmap=plt.cm.gray)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*get_frame_corners(shot, camera)), color='r')
    plt.xlim(extents[0:2])
    plt.ylim(extents[2:4])

    # H_alpha timeseries plot
    time_ha2, ha2 = get_time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    # Line average density plot
    time_dens, dens = time_crop(get_time_dens(shot), time)
    ax2 = plt.subplot(gs[2])#.locator_params(tight=True, nbins=6)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])

    # Slider and button settings
    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'Frame', 0, frame_count, valinit=0)
    slider.drawon = True
    slider.valfmt = '%d'
    playbutton_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    playbutton = Button(playbutton_area, 'Play')

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
                                       frames=frame_count, interval=0, blit=False)

    playbutton.on_clicked(play)
    plt.show()


if __name__ == '__main__':
    # Cziegler: 1101209014 with apd_array 
    shot = 1150611004 #1150528015   
    camera = 'phantom2' # try phantom
    time = get_gpi_series(shot, camera, 'time')
    frames = flip_horizontal(get_gpi_series(shot, camera, 'frames'))
    frames = subtract_average(frames, 5)

    animate_video(shot, camera, time, frames)
    #slide_frames(shot, camera, time, frames)

