from scipy.io.idl import readsav
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import acquire
import signal
import view


def mark_filament(base, overlay):
    newbase = np.ravel(base.copy())
    newoverlay = np.ravel(overlay)
    for i in range(len(newbase)):
        if newoverlay[i] > 0: newbase[i] = np.nan
    return newbase.reshape((64, 64))


sav = readsav('/usr/local/cmod/codes/spectroscopy/ir/FLIR/fl_images_1150611004.sav')
# Get field line projection images
fls = sav.fl_image
#gpi.slide_corr(fls, (0, 0))

#frames = gpi.flip_horizontal(gpi.get_gpi_series(1150611004, 'phantom2', 'frames'))
#frames = gpi.edge_filter(gpi.subtract_average(frames[:2048], 5))
#frames = gpi.kill_sobel_edges(frames)
subs = acquire.video(shot, camera, sub=5, sobel=True)
frames = subs[:2048]

gpi_index = 1028
xcorrs = np.zeros(fls.shape[0])
for i, fl in enumerate(fls):
    xcorrs[i] = signal.cross_correlation(frames[gpi_index], fl) 
indices = np.argsort(xcorrs)

fl_r = []; fl_z = []
noframes = 256
for j, frame in enumerate(frames[:noframes]):
    xcorrs = np.zeros(fls.shape[0])
    for i, fl in enumerate(fls):
        xcorrs[i] = signal.cross_correlation(frame, fl) 
    indices = np.argsort(xcorrs)
    r = sav.fieldline_r[indices[-1]]; z = sav.fieldline_z[indices[-1]]
    fl_r.append(r); fl_z.append(z)
    if not j % 10: print j, '/', noframes

view.plot_field_lines(fl_r, fl_z)

#plt.figure(); plt.plot(xcorrs); plt.show()

# Plot field line overlaid on GPI frame
fig, ax = plt.subplots()
cmap = plt.cm.gray; cmap.set_bad((1, 0, 0, 1))
im_over = plt.imshow(mark_filament(frames[gpi_index], fls[indices[-1]]), cmap=cmap, origin='bottom')

# Slider settings
fl_slide_area = plt.axes([0.17, 0.02, 0.65, 0.03])
fl_slider = Slider(fl_slide_area, 'Field line', 0, len(fls)-1, valinit=0)
fl_slider.valfmt = '%d'
gpi_slide_area = plt.axes([0.17, 0.06, 0.65, 0.03])
gpi_slider = Slider(gpi_slide_area, 'GPI frame', 0, len(frames)-1, valinit=0)
gpi_slider.valfmt = '%d'

def update_base(val):
    global gpi_index; global indices
    gpi_index = val
    for i, fl in enumerate(fls):
        xcorrs[i] = signal.cross_correlation(frames[val], fl) 
    indices = np.argsort(xcorrs)
    update_overlay(-1)

def update_overlay(val):
    global gpi_index; global indices
    im_over.set_array(mark_filament(frames[gpi_index], fls[indices[val]]))
    fig.canvas.draw_idle()

fl_slider.on_changed(update_overlay)
gpi_slider.on_changed(update_base)
plt.subplots_adjust(bottom=0.1)
plt.show()

