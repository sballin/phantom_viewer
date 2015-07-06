import matplotlib.pyplot as plt
import eqtools
import gpi

fl_r=[.547, .54, .555, .565] 
fl_z=[-.48, -.49, -.47, -.47] 

efit_tree = eqtools.CModEFIT.CModEFITTree(1150611004)
efit_times = efit_tree.getTimeBase()
rlcfs = efit_tree.getRLCFS()
zlcfs = efit_tree.getZLCFS()
machine_x, machine_y = efit_tree.getMachineCrossSectionFull()
corners = gpi.get_frame_corners(1150611004, 'phantom2')
corners_r, corners_z = [c[0] for c in corners], [c[1] for c in corners]

plt.figure()
plt.plot(fl_r, fl_z, 'bo')
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

