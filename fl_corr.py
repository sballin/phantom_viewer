import os
from scipy.io.idl import readsav
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import acquire
import signals
import view


def make_fls(shot):
    os.system("idl -e '.r /usr/local/cmod/codes/spectroscopy/ir/FLIR/create_periscope_view_Xpt.pro,'")


def overlay_filament(base, overlay):
    newbase = np.ravel(base.copy())
    newoverlay = np.ravel(overlay)
    for i in range(len(newbase)):
        if newoverlay[i] > 0: newbase[i] = np.nan
    return newbase.reshape((64, 64))


def get_best_fl(shot, fls, time):
    xr, xz = acquire.x_pt_r_z(shot)
    xrz = xrz[40] #######JUNK
    distances = np.zeros(fls.shape[0])
    for i, f in enumerate(fls):
        distances[i] = np.sqrt((sav.fieldline_r[i]-xrz[0])**2+(sav.fieldline_z[i]-xrz[1])**2)
    plt.figure()
    plt.imshow(fls[np.argmin(distances)], origin='bottom', cmap=plt.cm.gray)
    print sav.fieldline_r[np.argmin(distances)], sav.fieldline_z[np.argmin(distances)]
    plt.show()
 

def show_x_pt_fl(shot):
    xr, xz = acquire.x_pt_r_z(shot)
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


def plot_fl_slider(shot, sav):
    """
    Plot field line overlaid on GPI frame alterable with slider.
    """
    fls = sav.fl_image
    fl_r = sav.fieldline_r; fl_z = sav.fieldline_z
    subs = acquire.video(shot, 'phantom2', sub=5, sobel=False)
    frames = subs[:2048]
    gpi_index = 1028
    xcorrs = np.zeros(fls.shape[0])
    for i, fl in enumerate(fls):
        xcorrs[i] = signals.cross_correlation(frames[gpi_index], fl) 
    indices = np.argsort(xcorrs)

    fig, ax = plt.subplots()
    plt.subplot(121)
    cmap = plt.cm.gray; cmap.set_bad((1, 0, 0, 1))
    im_over = plt.imshow(overlay_filament(frames[gpi_index], fls[indices[-1]]), cmap=cmap, origin='bottom')

    rlcfs = acquire.rlcfs(shot)
    zlcfs = acquire.zlcfs(shot)
    machine_x, machine_y = acquire.machine_cross_section()

    ax1 = plt.subplot(122)
    plt.scatter(fl_r, fl_z, marker='o')
    plt.plot(machine_x, machine_y, color='gray')
    plt.plot(rlcfs[60], zlcfs[60])
    plt.axis('equal')
    plt.xlim([.47, .64]); plt.ylim([-.55, -.25])
    f, = plt.plot(fl_r[indices[-1]], fl_z[indices[-1]], 'ro')
    
    # Slider and button settings
    fl_slide_area = plt.axes([0.17, 0.02, 0.65, 0.03])
    fl_slider = Slider(fl_slide_area, 'Fit', 0, len(fls)-1, valinit=0)
    fl_slider.valfmt = '%d'
    gpi_slide_area = plt.axes([0.17, 0.06, 0.65, 0.03])
    gpi_slider = Slider(gpi_slide_area, 'GPI frame', 0, len(frames)-1, valinit=0)
    gpi_slider.valfmt = '%d'
    forward_button_area = plt.axes([0.95, 0.06, 0.04, 0.04])
    forward_button = Button(forward_button_area, '>')
    back_button_area = plt.axes([0.95, 0.01, 0.04, 0.04])
    back_button = Button(back_button_area, '<')
    
    def update_base(val):
        global gpi_index; global indices
        gpi_index = val
        for i, fl in enumerate(fls):
            xcorrs[i] = signals.cross_correlation(frames[val], fl) 
        indices = np.argsort(xcorrs)
        update_overlay(-1)
    
    def update_overlay(val):
        global gpi_index; global indices
        ax1.set_title('R: %.5f, Z: %.5f' % (sav.fieldline_r[indices[val]], sav.fieldline_z[indices[val]]), color='r')
        im_over.set_array(overlay_filament(frames[gpi_index], fls[indices[val]]))
        f.set_xdata(fl_r[indices[val]])
        f.set_ydata(fl_z[indices[val]])
        fig.canvas.draw_idle()

    def forward(event):
        global gpi_index
        gpi_index += 1
        update_base(gpi_index)

    def backward(event):
        global gpi_index
        gpi_index -= 1
        update_base(gpi_index)
    
    fl_slider.on_changed(update_overlay)
    gpi_slider.on_changed(update_base)
    forward_button.on_clicked(forward)
    back_button.on_clicked(backward)

    plt.subplots_adjust(bottom=0.1)
    plt.show()
        

if __name__ == '__main__':
    shot = 1150611004
    #show_x_pt_fl(shot)
    
    # Get field line projection images
    sav = readsav('/usr/local/cmod/codes/spectroscopy/ir/FLIR/fl_images_1150611004.sav')
    fls = sav.fl_image
    #fls_cross_section_plot(sav)
    plot_fl_slider(shot, sav)

