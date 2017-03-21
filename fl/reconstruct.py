import numpy as np
from phantom_viewer import acquire
from phantom_viewer import process
from phantom_viewer.fl import make_fl_images
from phantom_viewer.fl import view
import scipy.io.idl 
import scipy.sparse
import invert_fil_sart
import glob
import os
import matplotlib.pyplot as plt


def reconstruct_sparse(shot, fl_sav):
    """
    Reconstruct Phantom camera frames using given synthetic field line images.
    Provides speedup when dealing with large matrices. Source: James Harrison,
    12/2016.
    
    Args:
        shot: [int] shot number
        fl_sav: [scipy.io.idl.readsav object] containing field line images
    """
    x = 64
    y = 64
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    fl_images = fl_sav.fl_image
    fl_r = fl_sav.fieldline_r
    fl_z = fl_sav.fieldline_z
    geomat_rows = np.zeros(len(fl_images)*x*y, dtype=np.uint16)
    geomat_cols = np.zeros(len(fl_images)*x*y, dtype=np.uint16)
    geomat_vals = np.zeros(len(fl_images)*x*y, dtype=np.float16)
    geomat_shape = (x*y, len(fl_images))

    ctr = 0
    for i, fl in enumerate(fl_images):
        fl = fl.flatten()
        indices = np.where(fl/float(np.max(fl)) > 1.e-9)
        geomat_cols[ctr:ctr+len(indices[0])] = i

        geomat_rows[ctr:ctr+len(indices[0])] = indices[0]
        geomat_vals[ctr:ctr+len(indices[0])] = fl[indices].astype(np.float16)
        ctr += len(indices[0])
        indices = None
        
    geomat_rows = geomat_rows[0:ctr]
    geomat_cols = geomat_cols[0:ctr]
    geomat_vals = geomat_vals[0:ctr]
    
    a = scipy.sparse.csc_matrix(geomat_vals)
    b = scipy.sparse.csc_matrix(frames[0])
    invert_fil_sart.invert_sart(a, b, lam_start=1)
    
    
def test_reconstruct():
    """
    Reconstruct one frame using scipy NNLS. 
    """
    shot = 1150611004
    fl_images = np.load('../cache/fl_images_Xpt_1150611004_00.npy')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)

    geomatrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
    target = frames[6249].flatten()
    ems, rnorm = scipy.optimize.nnls(geomatrix, target)
    reconstruction = geomatrix.dot(ems).reshape((64,64))
    
    plt.figure()
    plt.subplot(121)
    plt.imshow(frames[6249],cmap=plt.cm.gist_heat, origin='bottom')
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(reconstruction,cmap=plt.cm.gist_heat, origin='bottom')
    plt.colorbar()
    plt.show()


def test_reconstruct_smooth():
    """
    Reconstruct one frame using scipy NNLS with least squares smoothing factor.
    """
    shot = 1150611004
    fl_images = np.load('../cache/fl_images_Xpt_1150611004_00.npy')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)

    lamda = 5000
    geomatrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
    geomatrix_smooth = np.concatenate((geomatrix, lamda*np.identity(len(fl_images))))
    target = frames[6249]
    target_smooth = np.concatenate((target.flatten(), np.zeros(len(fl_images))))
    return scipy.optimize.nnls(geomatrix_smooth, target_smooth)
    

def test_smoothing(shot, fl_sav):
    """
    Show results of using various smoothing parameters.
    """
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)

    for lamda in [0.1, 1, 10, 100, 1000, 5000, 10000]:
        view.slide_reconstruct(shot, fl_sav, smoothing_param=lamda)


def test_julia():
    shot = 1150611004
    smoothing_param = 100

    frames = acquire.video(shot, 'phantom2', sub=20)[:100]
    fl_images = np.load('../cache/fl_images_Xpt_{}_{:02d}.npy'.format(shot, 0))

    if smoothing_param:
        fl_matrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
        fl_matrix = np.concatenate((fl_matrix, smoothing_param*np.identity(len(fl_images))))
        frames_flattened = np.array([np.concatenate((f.flatten(), np.zeros(len(fl_images)))) for f in frames])
    else:
        fl_matrix = np.transpose(np.array([fl.flatten().astype(float) for fl in fl_images]))
        frames_flattened = np.array([f.flatten().astype(float) for f in frames])
    np.save('../cache/fl_matrix_Xpt_{}_{}.npy'.format(shot, 0), fl_matrix)
    np.save('../cache/frames_Xpt_{}_{}.npy'.format(shot, 0), frames_flattened)

    os.system('julia -p 8 nnls.jl {} {} {}'.format(shot, 1, smoothing_param))

    os.system('rm -rf ../cache/frames_Xpt_{}*.npy'.format(shot))
    os.system('rm -rf ../cache/fl_matrix_Xpt_{}*.npy'.format(shot))


def write_nnls_reconstruction(shot, smoothing_param=100):
    """
    Reconstruct entire shot and save field line emissivity profiles.

    Args:
        shot: [int] shot number
    """
    # Change working directory to write cache files to correct locations
    old_working_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Get all required data
    frames = acquire.video(shot, 'phantom2', sub=20)
    phantom_times = acquire.gpi_series(shot, 'phantom2', 'time')

    # Group frames by nearest EFIT timestep
    efit_times = [round(t, 3) for t in acquire.times_efit(shot)]
    phantom_efit_times = [efit_times[process.find_nearest(efit_times, t, ordered=True)] for t in phantom_times] 
    efit_times = sorted([round(t, 3) for t in set(phantom_efit_times)])
    frames_grouped = [[] for i in range(len(efit_times))]
    for i, frame in enumerate(frames):
        frames_grouped[efit_times.index(phantom_efit_times[i])].append(frame)

    # Load/calc+save field line images 
    fl_files = sorted(glob.glob('../cache/fl_images_Xpt_{}*'.format(shot)))
    if fl_files:
        print "STATUS: Loaded field line image files"
    else:
        make_fl_images.write(shot, efit_times)
        fl_files = sorted(glob.glob('../cache/fl_images_Xpt_{}*'.format(shot)))

    # Save field line images and phantom frames for julia
    #efit_times = [efit_times[0]] # TEMPORARY PLS REMOVE ME LATER
    for i, time in enumerate(efit_times):
        fl_images = np.load('../cache/fl_images_Xpt_{}_{:02d}.npy'.format(shot, i))
        if smoothing_param:
            fl_matrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
            fl_matrix = np.concatenate((fl_matrix, smoothing_param*np.identity(len(fl_images))))
            frames_flattened = np.array([np.concatenate((f.flatten(), np.zeros(len(fl_images)))) for f in frames_grouped[i]])
        else:
            fl_matrix = np.transpose(np.array([fl.flatten().astype(float) for fl in fl_images]))
            frames_flattened = np.array([f.flatten().astype(float) for f in frames_grouped[i]])
        np.save('../cache/fl_matrix_Xpt_{}_{}.npy'.format(shot, i), fl_matrix)
        np.save('../cache/frames_Xpt_{}_{}.npy'.format(shot, i), frames_flattened)

    # Run julia and wait until completion
    print "STATUS: Invoking Julia for NNLS reconstruction"
    os.system('julia -p 8 nnls.jl {} {} {}'.format(shot, len(efit_times), smoothing_param))
    # Delete files saved for julia
    os.system('rm -rf ../cache/frames_Xpt_{}*.npy'.format(shot))
    os.system('rm -rf ../cache/fl_matrix_Xpt_{}*.npy'.format(shot))

    # Do not remove without removing old os.chdir call
    os.chdir(old_working_dir)


def main():
    for shot in [1160505008, 1160505011]:
        print 'STATUS: Working on shot {}'.format(shot)
        write_nnls_reconstruction(shot)
        acquire.Database().purge()


if __name__ == '__main__':
    main()
