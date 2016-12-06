import sys
import os
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.gridspec as gridspec
import numpy as np
import eqtools
import acquire
import process


def animate_gpi(shot, camera='phantom2', sub=20, blur=0, sobel=False, 
                interval=20, skip=1):
    """
    Animate frames while displaying last closed flux surface and relevant plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        sub: number of frames to use in average subtraction
        sobel: whether to apply a Sobel filter
        interval: time between frames in ms
        skip: number of frames to skip, 1 means don't skip any
    """
    # Get time and frames
    time = acquire.gpi_series(shot, camera, 'time')
    frames = acquire.video(shot, camera, sub=sub, blur=blur, sobel=sobel)
    frame_count = frames.shape[0]

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = process.find_nearest(efit_times, time[0])
    gpi_extent = acquire.extent(shot, camera)

    # GPI, LCFS, and flux contour initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=gpi_extent, 
                    cmap=plt.cm.gray)
    im.get_axes().locator_params(nbins=6)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='r')
    plt.scatter(*zip(*acquire.frame_corners(shot, camera)), color='r')
    plt.xlim(gpi_extent[0:2])
    plt.ylim(gpi_extent[2:4])

    # H_alpha timeseries plot
    time_ha2, ha2 = acquire.time_ha2(shot)
    plt.subplot(gs[1]).locator_params(axis='y', nbins=4)
    plt.plot(time_ha2, ha2)
    plt.title('H-alpha')
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('[Units]')
    plt.xlim([time[0], time[-1]])
    plt.ylim([0, np.max([s for i, s in enumerate(ha2) 
                         if time[0] < time_ha2[i] < time[-1]])])

    # Line average density plot
    time_dens, dens = process.time_crop(acquire.time_dens(shot),time) 
    plt.subplot(gs[2]).locator_params(axis='y', nbins=4)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('[Units]')
    plt.xlim([time[0], time[-1]])

    def init():
        im.set_data(frames[0])
        efit_t_index = process.find_nearest(efit_times, time[0])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[0])
        vl2.set_xdata(time[0])
        return [im, l, vl1, vl2]
    
    def animate(i):
        i *= skip
        im.set_array(frames[i])
        efit_t_index = process.find_nearest(efit_times, i/float(frame_count)*(time[-1]-time[0]) + time[0])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, vl1, vl2]
    
    dim = frames.shape
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=dim[0], interval=interval, blit=True)
    plt.tight_layout(pad=1)
    plt.show()


def slide_gpi(shot, camera='phantom2', sub=20, blur=3, interval=0, 
              pixel_t_hist=None):
    """
    Slide through GPI frames while displaying last closed flux surface and 
    relevant timeseries plots.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplange GPI),
                'phantom2' (X-point GPI)
        sub: number of frames to use in average subtraction. 20 works well.
        blur: extent of Gaussian blur, 1, 2, or 3 advisable
        interval: delay in ms between frames during play
        pixel_t_hist: (x,y) to show time history of that pixel instead of 
        H-alpha.
    """
    time = acquire.gpi_series(shot, camera, 'time')
    frames = acquire.video(shot, camera)
    frame_count = frames.shape[0]
    phantom_extent = acquire.extent(shot, camera)
    efit_times = acquire.times_efit(shot)
    efit_t_index = process.find_nearest(efit_times, time[0])
    
    # Plot GPI and LCFS 
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25, hspace=.43, top=.96)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=phantom_extent, 
                    cmap=plt.cm.hot)
    plt.colorbar()
    im.get_axes().locator_params(nbins=6)
    plt.scatter(*zip(*acquire.frame_corners(shot, camera)), color='r')
    plt.xlim(phantom_extent[0:2])
    plt.ylim(phantom_extent[2:4])

    # Plot H_alpha or pixel timeseries
    plt.subplot(gs[1]).locator_params(axis='y', nbins=4)
    if pixel_t_hist:
        plt.title('Pixel time history')
        plt.plot(time, 
                 frames[:, pixel_t_hist[0], pixel_t_hist[1]].swapaxes(0, 1))
    else: 
        plt.title('H-alpha')
        time_ha2, ha2 = process.time_crop(acquire.time_ha2(shot), time)
        plt.plot(time_ha2, ha2)
    vl1 = plt.axvline(time[0], color='r')
    plt.ylabel('[Units]')
    plt.xlim([time[0], time[-1]])

    # Plot line average density 
    time_dens, dens = process.time_crop(acquire.time_dens(shot), time) 
    plt.subplot(gs[2]).locator_params(axis='y', nbins=4)
    plt.title('Line Average Density')
    plt.plot(time_dens, dens)
    vl2 = plt.axvline(time[0], color='r')
    plt.xlabel('Time (s)')
    plt.ylabel('[Units]')
    plt.xlim([time[0], time[-1]])

    curr_frame = 0

    def update(val):
        global frames, curr_frame
        val = int(val)
        curr_frame = val 
        curr_time = time[val]
        slider.valtext.set_text('%d (t=%.5f s)' % (val, curr_time))
        if val < frame_count: 
            im.set_array(frames[val])
            efit_t_index = process.find_nearest(efit_times, curr_time)
            vl1.set_xdata(curr_time)
            vl2.set_xdata(curr_time)
        fig.canvas.draw_idle()
 
    # Slider settings
    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'Frame', 0, frame_count, valinit=0)
    slider.drawon = True
    slider.valfmt = '%d'
    slider.on_changed(update)

    def init():
        global curr_frame
        curr_frame = int(curr_frame)
        start_frame = curr_frame
        im.set_data(frames[curr_frame])
        vl1.set_xdata(time[curr_frame])
        vl2.set_xdata(time[curr_frame])
        for i in xrange(curr_frame, curr_frame + 200, 1): 
            slider.set_val(i)
            vl1.set_xdata(time[i])
            vl2.set_xdata(time[i])
        slider.set_val(start_frame)
        return [im, l, vl1, vl2]

    def animate(i):
        im.set_array(frames[i])
        vl1.set_xdata(time[i])
        vl2.set_xdata(time[i])
        return [im, l, vl1, vl2]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=frame_count, interval=interval, 
                                       blit=False)

    def forward(event):
        global curr_frame
        curr_frame += 1
        slider.set_val(curr_frame)
        update(curr_frame)

    def backward(event):
        global curr_frame
        curr_frame -= 1
        slider.set_val(curr_frame)
        update(curr_frame)

    # Time button settings
    play_button_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    play_button = Button(play_button_area, 'Play')
    play_button.on_clicked(play)
    forward_button_area = plt.axes([0.56, 0.05, 0.04, 0.025])
    forward_button = Button(forward_button_area, '>')
    forward_button.on_clicked(forward)
    back_button_area = plt.axes([0.56, 0.015, 0.04, 0.025])
    back_button = Button(back_button_area, '<')
    back_button.on_clicked(backward)

    def filter(label):
        global frames
        if label == 'Orig': 
            frames = acquire.video(shot, camera)
        elif label == 'Sub %d' % sub:
            frames = acquire.video(shot, camera, sub=sub)
        elif label == 'Blur %d' % blur:
            frames = acquire.video(shot, camera, sub=sub, blur=blur)
        elif label == 'Sobel':
            frames = acquire.video(shot, camera, sub=sub, blur=blur, 
                                   sobel=True)
        update(curr_frame)
        im.autoscale()

    def recolor(event):
        im.autoscale()

    def cmap_change(label):
        global curr_frame
        if label == 'Red': im.cmap = plt.cm.hot
        if label == 'Gray': im.cmap = plt.cm.gray
        update(curr_frame)
        im.autoscale()

    # Image button settings
    left = .3
    bottom = .79
    filter_radio_area = plt.axes([left, bottom, .1, .12])
    filter_radio = RadioButtons(filter_radio_area, ('Orig', 'Sub %d' % sub, 
                                                    'Blur %d' % blur, 
                                                    'Sobel'))
    filter_radio.on_clicked(filter)
    cmap_radio_area = plt.axes([left, bottom-.08, .1, .07])
    cmap_radio = RadioButtons(cmap_radio_area, ('Red', 'Gray'))
    cmap_radio.on_clicked(cmap_change)
    recolor_button_area = plt.axes([left, bottom-.13, 0.1, 0.04])
    recolor_button = Button(recolor_button_area, 'Recolor')
    recolor_button.on_clicked(recolor)

    # Black magic to fix strange global variable error
    filter('Orig')

    plt.show()


def animate_frames(frames, cmap=plt.cm.gray, disp=True, save=False, interval=5):
    """
    Animate given frames.
    Parameters
        frames: NumPy array of frames with dimension (frame count, x pixels,
                y pixels)
        cmap: colormap e.g. plt.cm.hot, plt.cm.RdBu_r
        disp: whether to display animated figure
        save: whether to save animated figure to 'out.gif'
        interval: time interval between frames in ms
    """
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0)
    im = plt.imshow(frames[0], origin='lower', cmap=cmap)

    def init():
        im.set_data(frames[0])
        return [im]
    
    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    plt.axis('off')
    plt.tight_layout(pad=0, rect=(0, 0, 1, 1))
    pos = ax.get_position(); pos.y0 = 0.; ax.set_position(pos)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frames.shape[0], interval=interval, 
                                   blit=False)
    if disp: plt.show()
    if save: anim.save('out.gif', writer='imagemagick', fps=20)


def animate_shot(shot, camera='phantom2', cmap=plt.cm.gray, disp=True, 
                 save=False, interval=50):
    """
    Animate frames for the given shot and camera.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom2', 'raspi2'
        cmap: colormap e.g. plt.cm.hot, plt.cm.RdBu_r
        disp: whether to display animated figure
        save: whether to save animated figure to 'out.gif'
        interval: time interval between frames in ms
    """
    frames = acquire.video(shot, camera)
    animate_frames(frames, cmap=cmap, disp=disp, save=save, interval=interval)


def output_gif(frames):
    """ 
    Save a set of frames as out.gif by producing temporary images and putting 
    them together using ImageMagick.
    Parameters
        frames: NumPy array with dimension (frame count, x pixels, y pixels)
    """
    frames = frames[:, ::-1, :]
    for i, frame in enumerate(frames):
        plt.imsave('outframe_%06d.png' % i, frame, vmin=frame.min(), 
                   vmax=frame.max(), cmap=plt.cm.hot)
    os.system('convert outframe_*.png -layers optimize -delay 2/100 out.gif && rm outframe_*.png')


def output_vid_frames(frames):
    frames = frames[:, ::-1, :]
    os.system('mkdir out ; rm -f out/*')
    fig = plt.imshow(frames[0])
    fig.set_cmap('hot')
    plt.axis('off')
    for i, frame in enumerate(frames):
        fig.set_data(frame, vmin=frame.min(), vmax=frame.max())
        plt.savefig('out/image%05d.png' % i, bbox_inches='tight')

        #plt.imsave('out/image%05d.png' % i, frame, vmin=frame.min(), vmax=frame.max(), cmap=plt.cm.hot)


def output_frames(shot, start, end, traces=True):
    camera = 'phantom2'
    time = acquire.gpi_series(shot, camera, 'time')
    if type(start) is int:
        start_frame = start
        end_frame = end
    else:
        start_frame = process.find_nearest(time, start)
        end_frame = process.find_nearest(time, end)
    frames = acquire.video(shot, camera, sub=20) 

    folder_name = '%s-%d-%d' % (shot, start_frame, end_frame)
    os.system('mkdir %s' % folder_name)

    if not traces:
        plt.figure()
        im = plt.imshow(frames[start_frame], origin='lower', cmap=plt.cm.hot, 
                        vmin=0, vmax=130)
        #plt.axis('off')
        plt.colorbar()
        for i in xrange(start_frame, end_frame):
            im.set_data(frames[i])
            #im.autoscale()
            plt.savefig(folder_name + '/frame%05d.png' % i, bbox_inches='tight')
        return

    # LCFS data gathering
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    rlcfs = efit_tree.getRLCFS()
    zlcfs = efit_tree.getZLCFS()
    efit_t_index = process.find_nearest(efit_times, time[0])
    gpi_extent = acquire.extent(shot, camera)

    # GPI, LCFS initial plotting
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    #fig, ax = plt.subplots()
    plt.figure()
    plt.subplots_adjust(hspace=.5)
    #plt.subplots_adjust(bottom=0.25)
    plt.subplot(gs[0])
    im = plt.imshow(frames[0], origin='lower', extent=gpi_extent, 
                    cmap=plt.cm.hot)
    im.get_axes().locator_params(nbins=6)
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='c')
    plt.xlim(gpi_extent[0:2])
    plt.ylim(gpi_extent[2:4])

    time_ha2, ha2 = acquire.time_ha2(shot)
    ha_max = np.max([s for i, s in enumerate(ha2) 
                     if time[start_frame] < time_ha2[i] < time[end_frame]])
    ha_min = np.min([s for i, s in enumerate(ha2) 
                     if time[start_frame] < time_ha2[i] < time[end_frame]])
    time_dens, dens = acquire.time_dens(shot)
    dens_max = np.max(dens)
    dens_min = np.min(dens)

    for i in xrange(start_frame, end_frame):
        # H-alpha
        plt.subplot(gs[1]).locator_params(axis='y', nbins=4)
        plt.title('H-alpha')
        plt.plot(time_ha2, ha2)
        plt.xlim([time[start_frame], time[end_frame]])
        plt.ylim([ha_min, ha_max])
        v1 = plt.axvline(time[i], color='r')
        #plt.ylabel('[Units]')
        # Line average density
        plt.subplot(gs[2]).locator_params(axis='y', nbins=4)
        plt.title('Line Average Density')
        plt.plot(time_dens, dens)
        plt.xlim([.5, 1.5])
        plt.ylim([dens_min, dens_max])
        plt.xlabel('Time (s)')
        #plt.ylabel('[Units]')
        v2 = plt.axvline(time[i], color='r')

        im.set_data(frames[i])
        im.autoscale()
        plt.savefig(folder_name + '/frame%05d' % i)
        v1.remove()
        v2.remove()


def plot_field_lines(shot, fl_r, fl_z):
    """
    Plot given field lines along with LCFS and psi_cont surfaces for given shot.
    Parameters
        shot: int, shot number
        fl_r: array, field line r values
        fl_z: array, field line z values
    """
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    efit_times = efit_tree.getTimeBase()
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    efit_t_index = process.find_nearest(efit_times, time[0])

    rlcfs, zlcfs = acquire.lcfs_rz(shot)
    machine_x, machine_y = acquire.machine_cross_section()
    # corners = acquire.frame_corners(shot, 'phantom2')
    # corners_r, corners_z = [c[0] for c in corners], [c[1] for c in corners]
    time, psirz, psiext = acquire.time_flux_extent(shot)
    
    plt.figure()
    plt.plot(fl_r, fl_z, 'b.')
    plt.plot(plt.contour(psirz[efit_t_index], levels=np.arange(np.min(psirz), 
             np.max(psirz), .007), extent=psiext))
    plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], 'r-')
    plt.plot(machine_x, machine_y, color='gray')
    #plt.plot(1.020, -.265, 'go')
    #plt.plot(corners_r, corners_z, 'go')
    #plt.annotate('Aperture', (1.020, -.265))
    #plt.annotate('Field lines', (fl_r[-1], fl_z[-1]))
    #plt.annotate('View corners', corners[2])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.axis('equal')
    plt.show()
    

def slide_corr(frames, pixel, other_pixels=None):
    """
    Display correlation values for a pixel with a slider for lag time.
    Parameters
        frames: NumPy array of correlations with dimension (no. lags, x pixels,
                y pixels)
        pixel: (x, y) to mark with a circle
        other_pixels: array of (x, y) values of pixels to mark
    """
    frame_count = frames.shape[0]

    # Initial plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(bottom=0.25)
    im = plt.imshow(frames[frame_count/2], origin='bottom', 
                    cmap=plt.cm.RdBu_r, vmin=-1, vmax=1)
    circ = plt.Circle(pixel[::-1], radius=1, edgecolor='r', fill=False)
    ax.add_patch(circ)
    if other_pixels: 
        t_hists_r = [o[0] for o in other_pixels]
        t_hists_z = [o[1] for o in other_pixels]
    else:
        t_hists_r = acquire.gpi_series(1150611004, 'phantom2', 'hist_xpix')
        t_hists_z = acquire.gpi_series(1150611004, 'phantom2', 'hist_ypix')
    for pos in zip(t_hists_r, t_hists_z): 
        ax.add_patch(plt.Circle(pos, radius=1, edgecolor='b', fill=False))
    plt.colorbar()
    plt.title('Correlation for pixel (%d, %d)' % pixel)
    plt.xlabel('Pixel y coordinate')
    plt.ylabel('Pixel x coordinate')

    # Slider and button settings
    slide_area = plt.axes([0.10, 0.1, 0.65, 0.03])
    slider = Slider(slide_area, 'Lag', -frame_count/2, frame_count/2, valinit=0)
    slider.drawon = True
    slider.valfmt = '%d'
    play_button_area = plt.axes([0.45, 0.025, 0.1, 0.04])
    play_button = Button(play_button_area, 'Play')

    curr_frame = 0

    def update(val):
        global curr_frame
        curr_frame = val + frame_count/2
        slider.valtext.set_text('%d' % val)
        if curr_frame < frame_count: 
            im.set_array(frames[curr_frame])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    def init():
        global curr_frame
        curr_frame = int(curr_frame)
        im.set_data(frames[curr_frame])
        for i in xrange(curr_frame, curr_frame + 200, 2): 
            slider.set_val(i)
        return [im]

    def animate(i):
        im.set_array(frames[i])
        return [im]
    
    def play(event):
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=
                                       frame_count, interval=0, blit=False)

    play_button.on_clicked(play)
    plt.show()


if __name__ == '__main__':
    try:
        shot = int(sys.argv[1])
    except:
        shot = None

    if shot: 
        slide_gpi(shot)
    else:
        print 'Try help(slide_gpi) and help(animate_gpi) for inspiration.'
