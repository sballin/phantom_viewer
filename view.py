import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
import matplotlib.gridspec as gridspec
import numpy as np
import eqtools
import acquire
import process


def animate_video(shot, camera, sub=5, sobel=True):
    """
    Animate frames while displaying last closed flux surface and relevant plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        sub: number of frames to use in average subtraction
        frames: whether to apply a Sobel filter
    """
    # Get time and frames
    time = acquire.gpi_series(shot, camera, 'time')
    frames = acquire.video(shot, camera, sub=sub, sobel=sobel)
    frame_count = frames.shape[0]

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = process.find_nearest(efit_times, time[0])
    gpi_extent = acquire.extents(shot, camera)
    psirz, psiext = acquire.psi_contours(shot)

    # GPI, LCFS, and flux contour initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    ax0 = plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=gpi_extent, cmap=plt.cm.gray)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    c, = plt.plot(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz), np.max(psirz), .001), extent=psiext))
    print im, l, c
    plt.scatter(*zip(*acquire.frame_corners(shot, camera)), color='r')
    plt.xlim(gpi_extent[0:2])
    plt.ylim(gpi_extent[2:4])

    # H_alpha timeseries plot
    time_ha2, ha2 = acquire.time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    # Line average density plot
    time_dens, dens = process.time_crop(acquire.time_dens(shot), time)
    ax2 = plt.subplot(gs[2])#.locator_params(tight=True, nbins=6)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])

    def init():
        im.set_data(frames[0])
        efit_t_index = process.find_nearest(efit_times, time[0])
        #c.set_array(psirz[efit_t_index])
        #c.set_data(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz[0]), np.max(psirz[0]), .001), extent=psiext))
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[0])
        vl2.set_xdata(time[0])
        return [im, l, c, vl1, vl2]
    
    def animate(i):
        im.set_array(frames[i])
        efit_t_index = process.find_nearest(efit_times, i/float(frame_count)*(time[-1]-time[0]) + time[0])
        #c.set_array(psirz[efit_t_index])
        #c.set_data(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz[0]), np.max(psirz[0]), .001), extent=psiext))
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, c, vl1, vl2]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=50, blit=True)
    plt.tight_layout(pad=1)
    plt.show()


def slide_frames(shot, camera, sub=5, sobel=True):
    """
    Slide through GPI frames while displaying last closed flux surface and 
    relevant timeseries plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        sub: number of frames to use in average subtraction
        sobel: whether to apply a Sobel filter
    """
    time = acquire.gpi_series(shot, camera, 'time')
    frames = acquire.video(shot, camera, sub=sub, sobel=sobel)
    frame_count = frames.shape[0]

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = process.find_nearest(efit_times, time[0])
    extents = acquire.extents(shot, camera)
    
    # GPI and LCFS initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=extents, cmap=plt.cm.gray)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*acquire.frame_corners(shot, camera)), color='r')
    plt.xlim(extents[0:2])
    plt.ylim(extents[2:4])

    # H_alpha timeseries plot
    time_ha2, ha2 = acquire.time_ha2(shot)
    ax1 = plt.subplot(gs[1])#.locator_params(tight=True, nbins=6)
    plt.plot(time_ha2, ha2)
    plt.title('H alpha 2')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('Units?')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, 3])

    # Line average density plot
    time_dens, dens = process.time_crop(acquire.time_dens(shot), time)
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
            efit_t_index = process.find_nearest(efit_times, curr_time)
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


def animate_simple(shot, camera, cmap=plt.cm.gray, disp=True, save=False):
    """
    Animate the given frames.
    Parameters
        shot: int, shot number
        camera: string e.g. 'raspi2'
        cmap: colormap e.g. plt.cm.hot, plt.cm.RdBu_r
        disp: whether to display animated figure
        save: whether to save animated figure to 'out.gif'
    """
    frames = acquire.video(shot, camera)
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


def plot_field_lines(fl_r, fl_z):
    efit_tree = eqtools.CModEFIT.CModEFITTree(1150611004)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    machine_x, machine_y = efit_tree.getMachineCrossSectionFull()
    corners = gpi.get_frame_corners(1150611004, 'phantom2')
    corners_r, corners_z = [c[0] for c in corners], [c[1] for c in corners]
    
    plt.figure()
    plt.plot(fl_r, fl_z, 'b.')
    plt.plot(rlcfs[80], zlcfs[80], 'r-')
    plt.plot(machine_x, machine_y, color='gray')
    plt.plot(1.020, -.265, 'go')
    plt.plot(corners_r, corners_z, 'go')
    plt.annotate('Aperture', (1.020, -.265))
    plt.annotate('Field lines', (fl_r[-1], fl_z[-1]))
    plt.annotate('View corners', corners[2])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.axis('equal')
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
    t_hists_r = acquire.gpi_series(1150611004, 'phantom2', 'hist_xpix')
    t_hists_z = acquire.gpi_series(1150611004, 'phantom2', 'hist_ypix')
    for pos in zip(t_hists_r, t_hists_z): ax.add_patch(plt.Circle(pos, radius=1, edgecolor='b', fill=False))
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


if __name__ == '__main__':
    # Cziegler: 1101209014 with apd_array 
    shot = 1150528015   
    camera = 'phantom2' 
    animate_video(shot, camera, sub=0, sobel=False)
