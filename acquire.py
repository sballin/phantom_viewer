import sys
import MDSplus
import numpy as np
import process
import eqtools


class Database:
    vids = dict()

    def purge(self):
        for key in self.vids.keys(): del(self.vids[key])


def video(shot, camera, sub=0, blur=0, sobel=False): 
    key = str(shot) + camera
    try: out = Database().vids[key]
    except KeyError: 
        if camera == 'phantom' or camera == 'phantom2':
            out = Database().vids[key] \
                = process.flip_horizontal(gpi_series(shot, camera, 'frames'))
        else: out = Database().vids[key] = other_video(shot, camera)
    if sub: 
        key += 'sub%d' % sub
        try: out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] = process.subtract_min(out, sub)
    if blur:
        key += 'gauss%d' % blur
        try: out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] = process.gauss(np.copy(out), blur)
    if sobel: 
        key += 'sobel'
        try: out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] \
                = process.kill_sobel_edges(process.sobel(np.copy(out)))
    return out


def gpi_series(shot, camera, series_name):
    """
    Get specified GPI-related series.
    Parameters
        shot: int, shot number
        camera: string e.g. 'phantom' (outboard midplane GPI),
                'phantom2' (X-point GPI)
        series_name: string e.g. 'time' (timepoints), 'frames' (array of frames 
                     with dimension (frame count, y pixels, x pixels)
    Returns
        series: chosen with series_name
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    if series_name == 'time': 
        try: 
            series = tree.getNode('gpi.%s.t_hists' % camera).dim_of().data()
        except MDSplus._tdishr.TdiException: 
            print 'Time not loaded from t_hists'
            start = tree.getNode('gpi.%s.settings.trig_time' % camera).data()
            frame_rate = tree.getNode('gpi.%s.settings.frame_rate' 
                                      % camera).data()
            num_frames = tree.getNode('gpi.%s.settings.num_frames' 
                                      % camera).data()
            return np.arange(start, start + num_frames/float(frame_rate), 
                             1./frame_rate)
    else: 
        try: 
            series = tree.getNode('gpi.%s.%s' % (camera, series_name)).data()
        except MDSplus._tdishr.TdiException: 
            sys.exit('ERROR loading %s' % series_name)
    print 'Got %s for %s on shot %s from tree' % (series_name, camera, shot)
    return series 


def other_video(shot, camera='raspi2'):
    """
    Fetch frames from the given camera for the given shot.
    Parameters
        shot: int, shot number
        camera: string e.g. 'raspi2', 'matrox2', 'matrox3'
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    try: series = tree.getNode('video_daq.%s.frames' % camera).data()
    except MDSplus._tdishr.TdiException: 
        sys.exit('ERROR loading frames')
    return series 


def time_ha2(shot):
    """
    Get timepoints and H_alpha signal for specified shot.
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    node = tree.getNode('\\ha_2_bright')
    return node.dim_of().data(), node.data()


def time_dens(shot):
    """
    Get timepoints and line average density for specified shot.
    """
    tree = MDSplus.Tree('electrons', shot)
    node = tree.getNode('tci.results.nl_04')
    return node.dim_of().data(), np.array(node.data())/.6e20


def time_probe(shot):
    tree = MDSplus.Tree('edge', shot)
    node = tree.getNode('probes.fmp.id_01.p0.i_slow')
    return node.dim_of().data(), node.data()


def times_efit(shot):
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    return efit_tree.getTimeBase()


def frame_corners(shot, camera):
    """
    Return (R, Z) FOV corner locations for given camera in meters.
    """
    if camera == 'phantom2':
        return [[.5095, -.429], [.5825, -.431], [.583, -.357], [.51, -.355]]
    else:
        tree = MDSplus.Tree('spectroscopy', shot)
        # Convert R, Z coordinates from centimeters to meters
        return [np.array(tree.getNode('gpi.%s.image_pos.%s_corner' 
                                      % (camera, corner)).data())/100.
                for corner in ['br', 'tr', 'tl', 'bl']]


def extent(shot, camera):
    """
    Calculate R, Z border locations and return in matplotlib 'extent' order.
    """
    bl, br, tr, tl = tuple(frame_corners(shot, camera))
    rmin = np.mean((bl[0], tl[0]))
    rmax = np.mean((br[0], tr[0]))
    zmin = np.mean((bl[1], br[1]))
    zmax = np.mean((tl[1], tr[1]))
    return [rmin, rmax, zmin, zmax]


def x_pt_rz(shot):
    """
    Return (R, Z) coordinate tuples of the lower X-point in m.
    """
    tree = MDSplus.Tree('analysis', shot)
    rseps = np.array(tree.getNode('efit.results.a_eqdsk.rseps').data())/100.
    zseps = np.array(tree.getNode('efit.results.a_eqdsk.zseps').data())/100.
    return rseps[:, 0], zseps[:, 0]
                                   

def time_flux_extent(shot):
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    time = efit_tree.getTimeBase()
    flux = efit_tree.getFluxGrid()
    rgrid = efit_tree.getRGrid()
    zgrid = efit_tree.getZGrid()
    extent = [np.min(rgrid), np.max(rgrid), np.min(zgrid), np.max(zgrid)]
    return time, flux, extent


def lcfs_rz(shot):
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    return efit_tree.getRLCFS(), efit_tree.getZLCFS()


def machine_cross_section():
    efit_tree = eqtools.CModEFIT.CModEFITTree(1150611004)
    machine_x, machine_y = efit_tree.getMachineCrossSectionFull()
    return machine_x, machine_y
