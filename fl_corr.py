import os
from scipy.io.idl import readsav
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import acquire
import signals
import view


def make_fls(shot):
    os.system("idl -e '.r /usr/local/cmod/codes/spectroscopy/ir/FLIR/create_periscope_view_Xpt.pro,{shot}'".format(shot))


def overlay_filament(base, overlay):
    newbase = np.ravel(base.copy())
    newoverlay = np.ravel(overlay)
    for i in range(len(newbase)):
        if newoverlay[i] > 0: newbase[i] = np.nan
    return newbase.reshape((64, 64))


def plot_xpoint_and_fieldlines(shot, sav):
    """
    Show positions of X-point during given shot and field line R, Z positions
    for given .sav file.
    """
    xr, xz = acquire.x_pt_rz(shot)
    plt.figure()
    plt.scatter(sav.fieldline_r, sav.fieldline_z, color='b')
    plt.scatter(xr, xz, color='r')
    plt.show()


def fls_cross_section_plot(shot, sav):
    """
    Show best matching field lines on cross section plot.
    """
    frames = acquire.video(shot, 'phantom2', sub=5, sobel=False)
    fl_r = fl_z = []
    noframes = 128
    for j, frame in enumerate(frames[:noframes]):
        xcorrs = np.zeros(fls.shape[0])
        for i, fl in enumerate(fls):
            xcorrs[i] = signals.cross_correlation(frame, fl) 
        fl_r.append(sav.fieldline_r[np.argmax(xcorrs)])
        fl_z.append(sav.fieldline_z[np.argmax(xcorrs)])
        if not j % 10: print j, '/', noframes
    view.plot_field_lines(shot, fl_r, fl_z)


def make_colormap(seq):
    """
    Return a LinearSegmentedColormap
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


def plot_fl_slider(shot, sav):
    """
    Match field lines to GPI frames alterable with slider.
    """
    fls = sav.fl_image
    fl_r = sav.fieldline_r
    fl_z = sav.fieldline_z
    rlcfs, zlcfs = acquire.lcfs_rz(shot)
    machine_x, machine_y = acquire.machine_cross_section()
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    gpi_index = 0

    # Find cross-correlation scores between frames and field line images
    xcorrs = np.zeros(fls.shape[0])
    for i, fl in enumerate(fls):
        xcorrs[i] = signals.cross_correlation(frames[gpi_index], fl) 
    indices = np.argsort(xcorrs)

    # Interpolate field line cross-correlation scores over R, Z grid
    r_space = np.linspace(min(fl_r), max(fl_r), 100)
    z_space = np.linspace(min(fl_z), max(fl_z), 100)
    r_grid, z_grid = np.meshgrid(r_space, z_space)
    xcorr_grid = matplotlib.mlab.griddata(fl_r, fl_z, xcorrs, r_grid, z_grid, 
                                          interp='linear')

    # Plot camera image with field line overlay
    fig, ax = plt.subplots()
    plt.subplot(121)
    plasma_image = plt.imshow(frames[gpi_index], cmap=plt.cm.gray, 
                              origin='bottom')
    overlay_cmap = make_colormap([(1., 0., 0., 0.), (1., 0., 0., 1.)])
    fl_image = plt.imshow(fls[indices[-1]], cmap=overlay_cmap, origin='bottom',
                          alpha=0.8)
    plt.title('Phantom camera frame')

    # Plot field line R, Z data in context of machine
    ax1 = plt.subplot(122)
    xcorr_image = plt.pcolormesh(r_grid, z_grid, xcorr_grid)
    colorbar = plt.colorbar()
    colorbar.set_label('Cross-correlation')
    plt.plot(machine_x, machine_y, color='gray')
    plt.plot(rlcfs[60], zlcfs[60], color='fuchsia')
    plt.axis('equal')
    plt.xlim([.49, .62])
    plt.ylim([-.50, -.33])
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    f, = plt.plot(fl_r[indices[-1]], fl_z[indices[-1]], 'ro')
    
    # Slider and button settings
    fl_slide_area = plt.axes([0.17, 0.02, 0.65, 0.03])
    fl_slider = Slider(fl_slide_area, 'Cross-correlation', 0, len(fls)-1, 
                       valinit=0)
    fl_slider.valfmt = '%d'
    gpi_slide_area = plt.axes([0.17, 0.06, 0.65, 0.03])
    gpi_slider = Slider(gpi_slide_area, 'Phantom frame no.', 0, len(frames)-1,
                        valinit=0)
    gpi_slider.valfmt = '%d'
    forward_button_area = plt.axes([0.95, 0.06, 0.04, 0.04])
    forward_button = Button(forward_button_area, '>')
    back_button_area = plt.axes([0.95, 0.01, 0.04, 0.04])
    back_button = Button(back_button_area, '<')
    plt.tight_layout()

    def update_data(val):
        global gpi_index, indices, xcorr_grid
        gpi_index = int(val)
        for i, fl in enumerate(fls):
            xcorrs[i] = signals.cross_correlation(frames[gpi_index], fl) 
        xcorr_grid = matplotlib.mlab.griddata(fl_r, fl_z, xcorrs, r_grid, 
                                              z_grid, interp='linear')
        indices = np.argsort(xcorrs)
        update_images(-1)
    
    def update_images(val):
        global gpi_index, indices, xcorr_grid
        ax1.set_title('R: %.5f, Z: %.5f' % (sav.fieldline_r[indices[val]], 
                                            sav.fieldline_z[indices[val]]), 
                      color='r')
        plasma_image.set_array(frames[gpi_index])
        fl_image.set_array(fls[indices[val]])
        xcorr_image.set_array(xcorr_grid[:-1, :-1].ravel())
        f.set_xdata(fl_r[indices[val]])
        f.set_ydata(fl_z[indices[val]])
        fig.canvas.draw_idle()

    def forward(event):
        global gpi_index
        gpi_index += 1
        update_data(gpi_index)

    def backward(event):
        global gpi_index
        gpi_index -= 1
        update_data(gpi_index)
    
    fl_slider.on_changed(update_images)
    gpi_slider.on_changed(update_data)
    forward_button.on_clicked(forward)
    back_button.on_clicked(backward)

    plt.subplots_adjust(bottom=0.1)
    plt.show()
        

def main():
    shot = 1150611004 

    # Get field line projection images
    sav = readsav('fl_images/fl_images_1150611004_780ms_test6.sav')

    # Show field lines in context of machine
    #fls_cross_section_plot(shot, sav)

    # Match field lines to images
    plot_fl_slider(shot, sav)


if __name__ == '__main__':
    main()
