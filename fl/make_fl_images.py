import numpy as np
import pidly
import scipy.io.idl 
import scipy.ndimage.filters
import PIL.Image, PIL.ImageDraw
import os
import glob


def write_pxdata(shot, times):
    """
    Call IDL script for each relevant EFIT timestep in the given shot to write 
    field line pixel curves and starting R, Z data.

    Args:
        shot: [int] shot number
        time: [float array] EFIT times to use for reconstruction
    """
    print "STATUS: Invoking IDL to calculate field line projections"
    idl = pidly.IDL('/usr/local/cmod/bin/idlpy')
    idl('.r make_Xpt_basis_frames.pro')
    for i, time in enumerate(times):
        idl("make_Xpt_basis_frames,{},{},'{:02d}'".format(shot, time, i))
   

def make_image(xpix, ypix):
    """
    Return 64x64 image converting pixel values into smooth lines.

    Args:
        xpix: [float array] x pixel values
        ypix: [float array] y pixel values
    """
    # Set image parameters
    scale_factor = 3
    end_side = 64
    margin = 20

    # Prepare pixels for an image [scale_factor] times larger than the final
    # excluding (0,0) which always gets added by IDL
    pxpairs = zip(xpix, ypix)
    pxpairs = [(scale_factor * (p[0] + margin), scale_factor * (p[1] + margin))
               for p in pxpairs if not p[0] == p[1] == 0] 

    # Draw thick lines on large image, scale down, and apply Gaussian filter
    im = PIL.Image.new('L', (scale_factor * (end_side+2*margin), scale_factor * (end_side+2*margin)))
    draw = PIL.ImageDraw.Draw(im)
    draw.line(pxpairs, fill=128, width=4)
    del draw
    im = im.resize((end_side + 2*margin, end_side + 2*margin), PIL.Image.LANCZOS)
    image = np.asarray(im, dtype=np.uint8)
    image = scipy.ndimage.filters.gaussian_filter(image, 2)
    image = image[margin:-margin, margin:-margin]
    return image


def write(shot, times):
    """
    Write image files to disk.

    Args:
        shot: [int] shot number
        time: [float array] EFIT times to use for reconstruction
    """
    # Change working directory to write cache files to correct locations
    old_working_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load/calc+save field line trajectory data files
    fl_data_files = sorted(glob.glob('../cache/fl_data_Xpt_{}*'.format(shot)))
    if fl_data_files:
        print "STATUS: Loaded field line data files"
    else:
        print "STATUS: Writing field line data files"
        write_pxdata(shot, times)
        fl_data_files = sorted(glob.glob('../cache/fl_data_Xpt_{}*'.format(shot)))

    # Write field line image files to disk
    print "STATUS: Writing field line image files"
    for i, filename in enumerate(fl_data_files): 
        fl_data = scipy.io.idl.readsav(filename)
        images = [make_image(fl_data.fl_map.xpixfl[j], fl_data.fl_map.ypixfl[j])
                  for j in range(len(fl_data.fl_map.xpixfl))]

        np.save('../cache/fl_images_Xpt_{}_{:02d}.npy'.format(shot, i), images)

    # Do not remove without removing old os.chdir call
    os.chdir(old_working_dir)
