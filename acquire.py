import os
import numpy as np
import process
import eqtools
import MDSplus


class Database:
    vids = dict()

    def purge(self):
        for key in self.vids.keys(): del(self.vids[key])
        

def video(shot, camera='phantom2', sub=0, blur=0, sobel=False, cache=False): 
    key = '{}_{}'.format(shot, camera)
    
    # Try to get from memory cache
    try: 
        out = Database().vids[key]
    except KeyError: 
        # Try to get from disk cache
        try:
            code_directory = os.path.dirname(os.path.abspath(__file__))
            out = Database().vids[key] \
                = np.load('{}/cache/frames_{}.npy'.format(code_directory, key))
        except:
            # Download using MDSplus, mem-cache, optionally disk-cache, and serve
            if camera == 'phantom' or camera == 'phantom2':
                out = Database().vids[key] \
                    = process.flip_horizontal(gpi_series(shot, camera, 'frames'))
            else: 
                out = Database().vids[key] \
                    = other_video(shot, camera)
            if cache:
                code_directory = os.path.dirname(os.path.abspath(__file__))
                np.save('{}/cache/frames_{}.npy'.format(code_directory, key), out)
    
    # For processed videos, check whether already mem-cached, if not mem-cache and serve
    if sub: 
        key = '{}_sub{}'.format(key, sub)
        try: 
            out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] \
                = process.subtract_min(out, sub)
    if blur:
        key = '{}_gauss{}'.format(key, blur)
        try: 
            out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] \
                = process.gauss(np.copy(out), blur)
    if sobel: 
        key = '{}_sobel'.format(key)
        try: 
            out = Database().vids[key] 
        except KeyError: 
            out = Database().vids[key] \
                = process.kill_sobel_edges(process.sobel(np.copy(out)))
                
    return out
    
    
def get_mds(shot, tree='spectroscopy', node, dim_of=False, cache=True):
    code_directory = os.path.dirname(os.path.abspath(__file__))
    dim_of_flag = '_dimof_' if dim_of else ''
    cached_filename = '{}/cache/{}{}_{}.npy'.format(code_directory, node, dim_of_flag, shot)
    # Try to get from disk cache
    try:
        out = np.load(cached_filename)
    except:
        # Try to get from MDSplus server on-site
        try:
            tree = MDSplus.Tree('spectroscopy', shot)
            if dim_of:
                out = tree.getNode(node).dim_of().data()
            else:
                out = tree.getNode(node).data()
                
        # Try to get from MDSplus server offsite
        except:
            connection = MDSplus.Connection('localhost')
            
        if cache:
            np.save(cached_filename, out)
    return out
        

def gpi_series(shot, camera, series_name, cache=True):
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
    if series_name == 'time': 
        try: 
            series = get_mds(shot, 'gpi.%s.t_hists' % camera, dim_of=True, cache=cache)
        except MDSplus._tdishr.TdiException: 
            print 'Time not loaded from t_hists'
            start = get_mds(shot, 'gpi.%s.settings.trig_time' % camera, cache=cache)
            frame_rate = get_mds(shot, 'gpi.%s.settings.frame_rate' % camera, cache=cache)
            num_frames = get_mds(shot, 'gpi.%s.settings.num_frames' % camera, cache=cache)
            return np.arange(start, start + num_frames/float(frame_rate), 
                             1./frame_rate)
    else:
        series = get_mds(shot, 'gpi.%s.%s' % (camera, series_name), cache=cache)
    return series


def other_video(shot, camera='raspi2'):
    """
    Fetch frames from the given camera for the given shot.
    Parameters
        shot: int, shot number
        camera: string e.g. 'raspi2', 'matrox2', 'matrox3'
    """
    tree = MDSplus.Tree('spectroscopy', shot)
    return tree.getNode('video_daq.%s.frames' % camera).data()


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


def xpt_rz(shot):
    """
    Return (R, Z) coordinate tuples of the lower X-point in m.
    """
    tree = MDSplus.Tree('analysis', shot)
    rseps = np.array(tree.getNode('efit.results.a_eqdsk.rseps').data())/100.
    zseps = np.array(tree.getNode('efit.results.a_eqdsk.zseps').data())/100.
    return rseps[:, 0], zseps[:, 0]
                                   

def time_flux_extent(shot, highres=False):
    lowres_efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    if highres:
        efit_tree = MDSplus.Tree('efit19', shot)
        time = efit_tree.getNode('results.g_eqdsk.rbbbs').dim_of().data()
        flux = efit_tree.getNode('results.g_eqdsk.psirz').data()
    else:
        time = lowres_efit_tree.getTimeBase()
        flux = lowres_efit_tree.getFluxGrid()
    rgrid = lowres_efit_tree.getRGrid()
    zgrid = lowres_efit_tree.getZGrid()
    extent = [np.min(rgrid), np.max(rgrid), np.min(zgrid), np.max(zgrid)]
    return (time, flux, extent)


def lcfs_rz(shot, highres=False):
    if highres:
        efit_tree = MDSplus.Tree('efit19', shot)
        return (efit_tree.getNode('results.g_eqdsk.rbbbs').data(), 
                efit_tree.getNode('results.g_eqdsk.zbbbs').data())
    else:
        efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
        return efit_tree.getRLCFS(), efit_tree.getZLCFS()


def machine_cross_section():
    efit_tree = eqtools.CModEFIT.CModEFITTree(1150611004)
    machine_x, machine_y = efit_tree.getMachineCrossSectionFull()
    return machine_x, machine_y


def current_sign(shot):
    """
    Steve Wolfe's method taken from eqtools which didn't
    implement it for CMod.
    """
    efit_tree = eqtools.CModEFIT.CModEFITTree(shot)
    if np.mean(efit_tree.getIpMeas()) > 1e5: 
        return 1
    else:
        return -1
