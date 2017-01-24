"""
Created on Thu Dec 15 18:18:59 2016

@author: James
"""

import numpy as np
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from skimage.feature import blob_dog, blob_log, blob_doh

dat = np.load('out.npz')


sav = readsav('data/Xpt_fieldlines_1150611004_780ms.sav')
fl_images = sav.fl_image
fl_r = sav.fieldline_r
fl_z = sav.fieldline_z

dr = np.min(np.diff(fl_r)[np.diff(fl_r)>0])
dz = np.min(np.diff(fl_z)[np.diff(fl_z)>0])


# Generate a regular grid to interpolate the volues onto
meshr = np.linspace(np.min(fl_r),np.max(fl_r),(np.max(fl_r)-np.min(fl_r))/dr)
meshz = np.linspace(np.max(fl_z),np.min(fl_z),(np.max(fl_z)-np.min(fl_z))/dz)

mesh_r, mesh_z = np.meshgrid(meshr,meshz)


blobs_r = []
blobs_z = []
emi_total = mesh_r*0.0
for j in range(1000):
    arr = dat['arr_0'][j]

    emi = mesh_r*0.0
    # Map the data onto a regular 2D array by filling the corners with zeros, and
    # not using interpolation
    for i in np.arange(len(fl_r)):
        dist = (mesh_r-fl_r[i])**2+(mesh_z-fl_z[i])**2
        if np.min(dist) < 0.001:
            indx = np.unravel_index(np.argmin((mesh_r-fl_r[i])**2+(mesh_z-fl_z[i])**2),mesh_r.shape)
            emi[indx] = arr[i]
            
    # Reshape the array into a 2D image
    emi = np.reshape(emi,mesh_r.shape)

    # Run a blob detector on the image
    det = blob_log(emi,threshold=np.std(emi)*3)
    emi_total += emi
    if det.shape[0] > 0:
        r = meshr[det[:,1].astype(int)]
        z = meshz[(det[:,0]).astype(int)]
        blobs_r.extend(r)
        blobs_z.extend(z)
        print j

plt.figure()
plt.imshow(emi_total, interpolation='nearest', extent=(np.min(mesh_r),np.max(mesh_r),np.min(mesh_z),np.max(mesh_z)),cmap='Greys')
plt.scatter(blobs_r, blobs_z, color='r')
plt.show()

# # Plot the resulting image
# plt.figure()
# plt.imshow(emi,interpolation='nearest',extent=(np.min(mesh_r),np.max(mesh_r),np.min(mesh_z),np.max(mesh_z)))
# for i in np.arange(len(det)):
    # r = meshr[det[i][1]]
    # z = meshz[-det[i][0]-1]
#     plt.plot(r,z,'og')

# plt.colorbar()
# plt.xlabel('R (m)')
# plt.ylabel('Z (m)')
# plt.show()

# plt.figure()
# plt.imshow(emi,interpolation='nearest',cmap='Greys')
# for i in np.arange(len(det)):
#     plt.plot(det[i][1],det[i][0],'og')
