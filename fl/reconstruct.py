import numpy as np
from phantom_viewer import acquire
from phantom_viewer import process
from phantom_viewer.fl import make_fl_images
import scipy.io.idl 
import scipy.sparse
import invert_fil_sart
import glob
import os


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
    sav = scipy.io.idl.readsav('cache/Xpt_fieldlines_1150611004_780ms.sav')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    fl_images = sav.fl_image

    lamda = 5000
    geomatrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
    geomatrix_smooth = np.concatenate((geomatrix, lamda*np.identity(len(fl_images))))
    target = frames[6249]
    target_smooth = np.concatenate((target.flatten(), np.zeros(len(fl_images))))
    np.save('cache/matrix.npy', geomatrix_smooth)
    np.save('cache/target.npy', target_smooth)
    return scipy.optimize.nnls(geomatrix_smooth, target_smooth)
    

def test_smoothing(shot, fl_sav):
    """
    Show results of using various smoothing parameters.
    """
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)

    for lamda in [0.1, 1, 10, 100, 1000, 5000, 10000]:
        slide_reconstruct(shot, fl_sav, smoothing_param=lamda)


def write_nnls_reconstruction(shot):
    """
    Reconstruct entire shot and save field line emissivity profiles.

    Args:
        shot: [int] shot number
    """
    # Change working directory to write cache files to correct locations
    old_working_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Get all required data
    frames = acquire.video(shot, 'phantom2', sub=20, sobel=False)
    phantom_times = acquire.gpi_series(shot, 'phantom2', 'time')

    # Group frames by nearest EFIT timestep
    efit_times = [round(t, 3) for t in acquire.times_efit(shot)]
    phantom_efit_times = [efit_times[process.find_nearest(efit_times, t)] for t in phantom_times] 
    efit_times = sorted([round(t, 3) for t in set(phantom_efit_times)])
    frames_grouped = [[] for i in range(len(efit_times))]
    for i, frame in enumerate(frames):
        frames_grouped[efit_times.index(phantom_efit_times[i])].append(frame)

    # remove in production!!!
    #for i, group in enumerate(frames_grouped):
    #    frames_grouped[i] = group[::700]

    # Load/calc+save field line images 
    fl_files = sorted(glob.glob('../cache/fl_images_Xpt_{}*'.format(shot)))
    if fl_files:
        print "STATUS: Loaded field line image files"
    else:
        make_fl_images.write(shot, efit_times)
        fl_files = sorted(glob.glob('../cache/fl_images_Xpt_{}*'.format(shot)))

    # Save field line images and phantom frames for julia
    for i, time in enumerate(efit_times):
        fl_images = np.load('../cache/fl_images_Xpt_{}_{:02d}.npy'.format(shot, i))
        fl_matrix = np.transpose(np.array([fl.flatten() for fl in fl_images]))
        np.save('../cache/fl_matrix_Xpt_{}_{}.npy'.format(shot, i), fl_matrix)
        frames_flattened = np.array([f.flatten() for f in frames_grouped[i]])
        np.save('../cache/frames_Xpt_{}_{}.npy'.format(shot, i), frames_flattened)

    # Run julia and wait until completion
    print "STATUS: Invoking Julia for NNLS reconstruction"
    os.system('julia -p 8 nnls.jl {} {}'.format(shot, len(efit_times)))
    # Delete saved phantom frame files because they take up a lot of space
    os.system('rm -rf ../cache/frames_Xpt_{}*.npy'.format(shot))

    # Do not remove without removing old os.chdir call
    os.chdir(old_working_dir)


def main():
    write_nnls_reconstruction(1150611004)


if __name__ == '__main__':
    main()
