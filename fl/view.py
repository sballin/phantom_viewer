import os
import scipy.io.idl 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
from phantom_viewer import acquire
from phantom_viewer import signals
from phantom_viewer import process
import scipy.optimize
import scipy.interpolate
import glob
import types


def make_colormap(seq):
    """
    Return a LinearSegmentedColormap.
    
    Args:
        seq: a sequence of floats and RGBA-tuples. The floats should be 
        increasing and in the interval (0,1).
    """
    seq = [(None,) * 4, 0.0] + list(seq) + [1.0, (None,) * 4]
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1, a1= seq[i - 1]
            r2, g2, b2, a2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
            cdict['alpha'].append([item, a1, a2])
    return matplotlib.colors.LinearSegmentedColormap('CustomMap', cdict)


def slide_correlation(shot, fl_sav):
    """
    Slide through Phantom camera frames using given synthetic field line images 
    to create a reconstruction of the original using cross-correlation to find
    the best match.
    
    Args:
        shot: [int] shot number
        fl_sav: [scipy.io.idl.readsav object] containing field line images
    """
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    frame_index = 0
    fls = fl_sav.fl_image
    fl_r = fl_sav.fieldline_r
    fl_z = fl_sav.fieldline_z
    rlcfs, zlcfs = acquire.lcfs_rz(shot)
    efit_times, flux, flux_extent = acquire.time_flux_extent(shot)
    efit_t_index = process.find_nearest(efit_times, time[0], ordered=True)
    machine_x, machine_y = acquire.machine_cross_section()

    # Find cross-correlation scores between frames and field line images
    xcorrs = np.zeros(fls.shape[0])
    for i, fl in enumerate(fls):
        xcorrs[i] = signals.cross_correlation(frames[frame_index], fl) 
    indices = np.argsort(xcorrs)

    # Interpolate field line cross-correlation scores over R, Z grid
    r_space = np.linspace(min(fl_r), max(fl_r), 100)
    z_space = np.linspace(min(fl_z), max(fl_z), 100)
    r_grid, z_grid = np.meshgrid(r_space, z_space)
    xcorr_grid = matplotlib.mlab.griddata(fl_r, fl_z, xcorrs, r_grid, z_grid, 
                                          interp='linear')

    # Plot camera image with field line overlay
    fig, ax = plt.subplots()
    fig.suptitle('Shot {}'.format(shot))

    plt.subplot(121)
    plt.title('Divertor camera view')
    plasma_image = plt.imshow(frames[frame_index], cmap=plt.cm.gray, 
                              origin='bottom')
    overlay_cmap = make_colormap([(1., 0., 0., 0.), (1., 0., 0., 1.)])
    fl_image = plt.imshow(fls[indices[-1]], cmap=overlay_cmap, origin='bottom',
                          alpha=0.8)
    plt.axis('off')

    # Plot field line R, Z data in context of machine
    plt.subplot(122)
    plt.title('Toroidal cross section')
    xcorr_image = plt.pcolormesh(r_grid, z_grid, xcorr_grid)
    colorbar = plt.colorbar()
    colorbar.set_label('Cross-correlation')
    plt.plot(machine_x, machine_y, color='gray')
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='fuchsia')
    plt.axis('equal')
    plt.xlim([.49, .62])
    plt.ylim([-.50, -.33])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    f, = plt.plot(fl_r[indices[-1]], fl_z[indices[-1]], 'ro')
    plt.contour(flux[efit_t_index], 100, extent=flux_extent)
    
    # Slider and button settings
    fl_slide_area = plt.axes([0.20, 0.02, 0.60, 0.03])
    fl_slider = Slider(fl_slide_area, 'Correlation rank', 0, len(fls)-1, 
                       valinit=0)
    fl_slider.valfmt = '%d'
    phantom_slide_area = plt.axes([0.20, 0.06, 0.60, 0.03])
    phantom_slider = Slider(phantom_slide_area, 'Camera frame', 0, len(frames)-1,
                        valinit=0)
    phantom_slider.valfmt = '%d'
    forward_button_area = plt.axes([0.95, 0.06, 0.04, 0.04])
    forward_button = Button(forward_button_area, '>')
    back_button_area = plt.axes([0.95, 0.01, 0.04, 0.04])
    back_button = Button(back_button_area, '<')

    def update_data(val):
        global frame_index, indices, xcorr_grid
        frame_index = int(val)
        for i, fl in enumerate(fls):
            xcorrs[i] = signals.cross_correlation(frames[frame_index], fl) 
        xcorr_grid = matplotlib.mlab.griddata(fl_r, fl_z, xcorrs, r_grid, 
                                              z_grid, interp='linear')
        indices = np.argsort(xcorrs)[::-1]
        update_images(0)
    
    def update_images(val):
        global frame_index, indices, xcorr_grid
        val = int(val)
        efit_t_index = process.find_nearest(efit_times, time[frame_index], ordered=True)
        plasma_image.set_array(frames[frame_index])
        plasma_image.autoscale()
        fl_image.set_array(fls[indices[val]])
        fl_image.autoscale()
        f.set_xdata(fl_r[indices[val]])
        f.set_ydata(fl_z[indices[val]])
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        xcorr_image.set_array(xcorr_grid[:-1, :-1].ravel())
        fig.canvas.draw_idle()

    def forward(event):
        global frame_index
        frame_index += 1
        update_data(frame_index)

    def backward(event):
        global frame_index
        frame_index -= 1
        update_data(frame_index)
    
    fl_slider.on_changed(update_images)
    phantom_slider.on_changed(update_data)
    forward_button.on_clicked(forward)
    back_button.on_clicked(backward)

    plt.tight_layout(rect=(0, .1, 1, .9))
    plt.show()
        
        
def slide_reconstruct(shot, fl_sav, smoothing_param=100):
    """
    Slide through Phantom camera frames using given synthetic field line images 
    to create a reconstruction of the original using non-negative least squares
    (NNLS) fitting. The NNLS method used here is slow, so there are long delays 
    between clicks and interface updates.
    
    Args:
        shot: [int] shot number
        fl_sav: [scipy.io.idl.readsav object] containing field line images
        smoothing_param: [float] least-squares smoothing parameter
    """
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    frame_index = 0
    fl_images = fl_sav.fl_image
    fl_r = fl_sav.fieldline_r
    fl_z = fl_sav.fieldline_z
    rlcfs, zlcfs = acquire.lcfs_rz(shot)
    efit_times, flux, flux_extent = acquire.time_flux_extent(shot)
    efit_t_index = process.find_nearest(efit_times, time[0], ordered=True)
    machine_x, machine_y = acquire.machine_cross_section()

    inversion_func = scipy.optimize.nnls
    geomatrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
    geomatrix_smooth = np.concatenate((geomatrix, 
                                       smoothing_param*np.identity(len(fl_images))))
    target = frames[0]
    target_smooth = np.concatenate((target.flatten(), np.zeros(len(fl_images))))
    inv = inversion_func(geomatrix_smooth, target_smooth)

    r_space = np.linspace(min(fl_r), max(fl_r), 100)
    z_space = np.linspace(min(fl_z), max(fl_z), 100)
    r_grid, z_grid = np.meshgrid(r_space, z_space)
    emissivity_grid = matplotlib.mlab.griddata(fl_r, fl_z, inv[0], r_grid,
                                               z_grid, interp='linear') 

    reconstructed = geomatrix.dot(inv[0])
    reconstructed = reconstructed.reshape((64,64))
    target = target.reshape((64,64))
    fig, ax = plt.subplots()
    
    plt.subplot(221)
    plt.title('Divertor camera view')
    plasma_image = plt.imshow(target, cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    plt.colorbar()
    
    plt.subplot(222)
    plt.title('Reconstruction')
    reconstruction_image = plt.imshow(reconstructed, cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    plt.colorbar()
    
    plt.subplot(223)
    plt.title('Reconstruction minus original')
    error_image = plt.imshow(target - reconstructed, cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    plt.colorbar()
    
    plt.subplot(224)
    plt.title('Toroidal cross section')
    emissivity_image = plt.pcolormesh(r_grid, z_grid, emissivity_grid)
    colorbar = plt.colorbar()
    colorbar.set_label('Relative emissivity')
    plt.axis('equal')
    plt.plot(machine_x, machine_y, color='gray')
    l, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='fuchsia')
    plt.xlim([.49, .62])
    plt.ylim([-.50, -.33])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.contour(flux[efit_t_index], 100, extent=flux_extent)

    phantom_slide_area = plt.axes([0.20, 0.02, 0.60, 0.03])
    phantom_slider = Slider(phantom_slide_area, 'Camera frame', 0, len(frames)-1,
                        valinit=0)
    phantom_slider.valfmt = '%d'
    forward_button_area = plt.axes([0.95, 0.06, 0.04, 0.04])
    forward_button = Button(forward_button_area, '>')
    back_button_area = plt.axes([0.95, 0.01, 0.04, 0.04])
    back_button = Button(back_button_area, '<')

    def update_data(val):
        global frame_index, emissivity_grid, reconstructed
        frame_index = int(val)
        target_smooth = np.concatenate((frames[frame_index].flatten(), np.zeros(len(fl_images))))
        inv = inversion_func(geomatrix_smooth, target_smooth)
        # inv = inversion_func(geomatrix, target)
        reconstructed = geomatrix.dot(inv[0]).reshape((64,64))
        emissivity_grid = matplotlib.mlab.griddata(fl_r, fl_z, inv[0], r_grid, 
                                              z_grid, interp='linear')
        plasma_image.set_array(frames[frame_index])
        reconstruction_image.set_array(reconstructed)
        emissivity_image.set_array(emissivity_grid[:-1, :-1].ravel())
        error_image.set_array(frames[frame_index]-reconstructed)
        error_image.autoscale()
        plasma_image.autoscale()
        reconstruction_image.autoscale()
        emissivity_image.autoscale()
        efit_t_index = process.find_nearest(efit_times, time[frame_index], ordered=True)
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        fig.canvas.draw_idle()

    def forward(event):
        global frame_index
        frame_index += 1
        update_data(frame_index)

    def backward(event):
        global frame_index
        frame_index -= 1
        update_data(frame_index)

    phantom_slider.on_changed(update_data)
    forward_button.on_clicked(forward)
    back_button.on_clicked(backward)

    plt.tight_layout(rect=(0, .1, 1, .9))
    plt.show()


def make_grid(rs, zs, r_grid, z_grid, values):
    try:
        return matplotlib.mlab.griddata(rs, zs, values, r_grid, z_grid, interp='linear') 
    except:
        points = np.array(zip(rs, zs), dtype=np.float)
        return scipy.interpolate.griddata(points, values, (r_grid, z_grid), method='linear', fill_value=0)


def cutoff_array(values, cutoff):
    if not cutoff:
        return values
    else:
        maxval = np.max(values)
        minval = cutoff/100.*maxval
        temp = np.clip(values-minval, 0, maxval)
        temp[np.nonzero(temp)] += minval
        return temp


def slide_reconstruction(shot, smoothing_param=100, save=False):
    """
    Slide through Phantom camera frames and their reconstructions from
    synthetic field line images calculated externally with Julia.
    
    Args:
        shot: [int] shot number
        fl_sav: [scipy.io.idl.readsav object] containing field line images
    """
    # Read cache files broken up by EFIT time segment
    old_working_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    emissivity_files = sorted(glob.glob('../cache/fl_emissivities_Xpt_{}_sp{}_*'.format(shot, smoothing_param)))
    fl_data_files = sorted(glob.glob('../cache/fl_data_Xpt_{}_*'.format(shot)))
    fl_image_files = sorted(glob.glob('../cache/fl_images_Xpt_{}_*'.format(shot)))
    if len(emissivity_files) == 0 or len(fl_data_files) == 0 or len(fl_image_files) == 0:
        raise Exception('Could not find files necessary to view reconstruction')
    emissivities = [np.load(f) for f in emissivity_files]
    fl_data = [scipy.io.idl.readsav(f) for f in fl_data_files]
    fl_images = [np.load(f) for f in fl_image_files]
    os.chdir(old_working_dir)

    frame_index = 0
    cutoff = 0
    times = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    fl_rs = [f.fieldline_r for f in fl_data]
    fl_zs = [f.fieldline_z for f in fl_data]
    rlcfs, zlcfs = acquire.lcfs_rz(shot)
    efit_times, flux, flux_extent = acquire.time_flux_extent(shot)
    machine_x, machine_y = acquire.machine_cross_section()
    efit_t_index = process.find_nearest(efit_times, times[0], ordered=True)
    efit_phantom_start_index = efit_t_index
    previous_segments_frame_count = [sum([len(e) for e in emissivities[:t_index]]) for t_index in range(len(emissivities))]

    geomatrices = [np.transpose(np.array([fl.flatten() for fl in fl_image_set])) for fl_image_set in fl_images]
    reconstructed = geomatrices[0].dot(emissivities[0][0])
    reconstructed = reconstructed.reshape((64,64))

    fl_r_all = np.concatenate(fl_rs)
    fl_z_all = np.concatenate(fl_zs)
    r_space = np.linspace(np.min(fl_r_all), np.max(fl_r_all), 100)
    z_space = np.linspace(np.min(fl_z_all), np.max(fl_z_all), 100)
    r_grid, z_grid = np.meshgrid(r_space, z_space)
    emissivity_grid = make_grid(fl_rs[0], fl_zs[0], r_grid, z_grid, emissivities[0][0])
    
    # Draw figure using first frame
    fig, ax = plt.subplots()
    title = plt.suptitle('Shot {} frame {}'.format(shot, 0))
    plt.tight_layout(rect=(0, .1, 1, .9))

    plt.subplot(221)
    plt.title('Divertor camera view')
    plasma_image = plt.imshow(frames[0], cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    #plt.colorbar()
    
    plt.subplot(222)
    plt.title('Reconstruction')
    reconstruction_image = plt.imshow(reconstructed, cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    #plt.colorbar()
    
    plt.subplot(223)
    plt.title('Original minus reconstruction')
    error_image = plt.imshow(frames[0] - reconstructed, cmap=plt.cm.gist_heat, origin='bottom')
    plt.axis('off')
    #plt.colorbar()
    
    plt.subplot(224)
    plt.title('Toroidal cross section')
    emissivity_image = plt.pcolormesh(r_grid, z_grid, emissivity_grid, cmap=plt.cm.plasma)
    #colorbar = plt.colorbar()
    #colorbar.set_label('Relative emissivity')
    plt.axis('equal')
    plt.plot(machine_x, machine_y, color='gray')
    l, = plt.plot(rlcfs[efit_phantom_start_index], 
                  zlcfs[efit_phantom_start_index], color='orange', alpha=0.5)
    plt.xlim([.49, .62])
    plt.ylim([-.50, -.33])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.contour(flux[efit_t_index], 300, extent=flux_extent, alpha=0.5)
    plt.gca().set_axis_bgcolor('black')

    if not save:
        phantom_slide_area = plt.axes([0.20, 0.02, 0.60, 0.03])
        phantom_slider = Slider(phantom_slide_area, 'Camera frame', 0, len(frames)-1,
                            valinit=0)
        phantom_slider.valfmt = '%d'
        cutoff_slider_area = plt.axes([0.20, 0.05, 0.60, 0.03])
        cutoff_slider = Slider(cutoff_slider_area, 'Cutoff', 0, 100, valinit=0)
        cutoff_slider.valfmt = '%d%%'
        forward_button_area = plt.axes([0.95, 0.06, 0.04, 0.04])
        forward_button = Button(forward_button_area, '>')
        back_button_area = plt.axes([0.95, 0.01, 0.04, 0.04])
        back_button = Button(back_button_area, '<')

    def update_frame(val):
        global efit_t_index, frame_index
        frame_index = int(val)
        cutoff = cutoff_slider.val

        # Recalculate things
        efit_t_index = process.find_nearest(efit_times, times[frame_index], ordered=True)
        efit_rel_index = efit_t_index - efit_phantom_start_index
        frame_rel_index = frame_index - previous_segments_frame_count[efit_rel_index]
        emissivity_cutoff = cutoff_array(emissivities[efit_rel_index][frame_rel_index], cutoff)
        emissivity_grid = make_grid(fl_rs[efit_rel_index], fl_zs[efit_rel_index], r_grid, z_grid, emissivity_cutoff) 
        reconstructed = geomatrices[efit_rel_index].dot(emissivity_cutoff).reshape((64,64))

        # Update canvas
        plasma_image.set_array(frames[frame_index])
        plasma_image.autoscale()
        reconstruction_image.set_array(reconstructed)
        reconstruction_image.autoscale()
        emissivity_image.set_array(emissivity_grid[:-1, :-1].ravel())
        emissivity_image.autoscale()
        error_image.set_array(frames[frame_index] - reconstructed)
        error_image.autoscale()
        l.set_xdata(rlcfs[efit_t_index])
        l.set_ydata(zlcfs[efit_t_index])
        title.set_text('Shot {} frame {}'.format(shot, frame_index))

    def forward(event):
        global frame_index
        phantom_slider.set_val(frame_index + 1)

    def backward(event):
        global frame_index
        phantom_slider.set_val(frame_index - 1)

    def update_cutoff(val):
        global frame_index
        update_frame(frame_index)

    if save:
        FFMpegWriter = animation.writers['ffmpeg']
        writer = FFMpegWriter(fps=15, bitrate=5000)
        anim = animation.FuncAnimation(fig, update_frame, frames=1000, interval=5, blit=False)
        anim.save('{}_sp{}.avi'.format(shot, smoothing_param), writer=writer)
        plt.close(fig)
    else:
        phantom_slider.on_changed(update_frame)
        cutoff_slider.on_changed(update_cutoff)
        forward_button.on_clicked(forward)
        back_button.on_clicked(backward)
        plt.show()


def animate_emissivity(shot, num_frames=1000, smoothing_param=100, highres=False):
        
    def setvisible(self, vis):
        for c in self.collections: c.set_visible(vis)
        
    def setanimated(self, ani):
        for c in self.collections: c.set_animated(ani)
        
    # Read cache files broken up by EFIT time segment
    old_working_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    emissivity_files = sorted(glob.glob('../cache/fl_emissivities_Xpt_{}_sp{}_*'.format(shot, smoothing_param)))
    emissivities = [np.load(f) for f in emissivity_files]
    fl_data_files = sorted(glob.glob('../cache/fl_data_Xpt_{}_*'.format(shot)))
    fl_data = [scipy.io.idl.readsav(f) for f in fl_data_files]
    os.chdir(old_working_dir)
    fl_rs = [f.fieldline_r for f in fl_data]
    fl_zs = [f.fieldline_z for f in fl_data]
    rlcfs, zlcfs = acquire.lcfs_rz(shot, highres=highres)

    # Get standard data
    times = acquire.gpi_series(shot, 'phantom2', 'time')
    machine_x, machine_y = acquire.machine_cross_section()
    if not num_frames:
        num_frames = len(times)
    
    # Set up indices and get flux data
    efit_times, flux, flux_extent = acquire.time_flux_extent(shot, highres=False)
    if highres:
        efit_times_highres, flux, _ = acquire.time_flux_extent(shot, highres=True)
        efit_t_index_highres = process.find_nearest(efit_times_highres, times[0], ordered=True)
    efit_t_index = process.find_nearest(efit_times, times[0], ordered=True)
    efit_phantom_start_index = efit_t_index
    previous_segments_frame_count = [sum([len(e) for e in emissivities[:t_index]]) for t_index in range(len(emissivities))]
    frame_index = 0

    # Set up grid on which to display emissivity profile
    fl_r_all = np.concatenate(fl_rs)
    fl_z_all = np.concatenate(fl_zs)
    r_space = np.linspace(np.min(fl_r_all), np.max(fl_r_all), 100)
    z_space = np.linspace(np.min(fl_z_all), np.max(fl_z_all), 100)
    r_grid, z_grid = np.meshgrid(r_space, z_space)
    emissivity_grid = matplotlib.mlab.griddata(fl_rs[0], fl_zs[0], emissivities[0][0], r_grid, z_grid, interp='linear') 

    # Plotting setup
    fig, ax = plt.subplots()
    emissivity_image = plt.pcolormesh(r_grid, z_grid, emissivity_grid, cmap=plt.cm.plasma)
    plt.axis('equal')
    plt.xlim([.49, .62])
    plt.ylim([-.50, -.33])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.gca().set_axis_bgcolor('black')
    machine, = plt.plot(machine_x, machine_y, color='gray')
    
    previous_efit_t_index = -1
    plots = [0]*num_frames
    
    # Compute plot elements for each frame
    for frame_index in range(num_frames):
        efit_t_index = process.find_nearest(efit_times, times[frame_index], ordered=True)
        efit_rel_index = efit_t_index - efit_phantom_start_index
        frame_rel_index = frame_index - previous_segments_frame_count[efit_rel_index]
        
        # Plot emissivity grid
        emissivity_grid = matplotlib.mlab.griddata(fl_rs[efit_rel_index], fl_zs[efit_rel_index], emissivities[efit_rel_index][frame_rel_index], r_grid, z_grid, interp='linear') 
        emissivity_image = plt.pcolormesh(r_grid, z_grid, emissivity_grid, cmap=plt.cm.plasma)
        emissivity_image.autoscale()
        
        # Plot LCFS and flux contours
        if highres:
            efit_t_index_highres = process.find_nearest(efit_times_highres, times[frame_index], ordered=True)
            if efit_t_index_highres > previous_efit_t_index:
                lcfs, = plt.plot(rlcfs[:, efit_t_index_highres], zlcfs[:, efit_t_index_highres], color='orange', alpha=0.5)
                flux_surfaces = plt.contour(flux[efit_t_index_highres], 300, extent=flux_extent, alpha=0.5)               
                # Monkey patch to enable animation for QuadContourSet objects
                flux_surfaces.set_visible = types.MethodType(setvisible, flux_surfaces)
                flux_surfaces.set_animated = types.MethodType(setanimated, flux_surfaces)
                flux_surfaces.axes = plt.gca()
                previous_efit_t_index = efit_t_index_highres
        elif efit_t_index > previous_efit_t_index:
            lcfs, = plt.plot(rlcfs[efit_t_index], zlcfs[efit_t_index], color='orange', alpha=0.5)
            flux_surfaces = plt.contour(flux[efit_t_index], 300, extent=flux_extent, alpha=0.5)       
            # Monkey patch to enable animation for QuadContourSet objects
            flux_surfaces.set_visible = types.MethodType(setvisible,flux_surfaces)
            flux_surfaces.set_animated = types.MethodType(setanimated,flux_surfaces)
            flux_surfaces.axes = plt.gca()
            previous_efit_t_index = efit_t_index
            
        text = 'Shot {} frame {}'.format(shot, frame_index)
        axtext = ax.text(90, 90, text)
        annotation = ax.annotate(text, xy=(0.35, 1.05), xycoords='axes fraction')
        plots[frame_index] = [lcfs, flux_surfaces, emissivity_image, machine, axtext, annotation]

    # Save animation to file
    FFMpegWriter = animation.writers['ffmpeg']
    writer = FFMpegWriter(fps=15, bitrate=5000)
    anim = animation.ArtistAnimation(fig, plots, interval=5, blit=False)
    anim.save('{}_sp{}_emissivity.avi'.format(shot, smoothing_param), writer=writer)
    plt.close(fig)

   
def make_videos():
    import gc
    import sys

    allshots = [1150611004, 1150717011, 1150625030, 1150820011, 1150929013, 1150929016, 
                1150923009, 1150923010, 1150923012, 1150923013, 1150923017, 1160505008,
                1160505011, 1150505013, 1150505014, 1150505015, 1150505016, 1150505017,
                1150505018, 1150505022, 1150505023, 1150505030]

    for shot in [1150625030, 1150820011, 1150923009, 1150923010, 1150923013, 1150923017, 1150929013, 1150929016]:
        try:
            animate_emissivity(shot, num_frames=100, highres=True)
        except:
            try:
                animate_emissivity(shot, num_frames=100, highres=False)
            except Exception as e: 
                print shot
                print str(e)
        gc.collect()


def main():
    shot = 1150625030
    slide_reconstruction(shot)


if __name__ == '__main__':
    main()
