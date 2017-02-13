#!/usr/bin/env python

"""
Example script running a field line tracer

Nick Walkden, May 2015
"""

from cyFieldlineTracer.cyfieldlineTracer import RK4Tracer, Euler2Tracer
import numpy as np
from scipy.interpolate import CubicSpline


def make_fls(shot, time):
    # Pick tracer, Euler2Tracer or RK4Tracer. tracer.eq shows equilibrium info.
    trace_method = RK4Tracer
    tracer = trace_method(machine='CMod', shot=shot, time=time, interp='quintic')

    # Set up field line grid and arrays to store traced paths
    rgrid = np.linspace(0.495, 0.495+0.115, 34)
    zgrid = np.linspace(-0.49, -0.49+0.14, 40)
    flr = [0]*(rgrid.shape[0]*zgrid.shape[0])
    flz = [0]*(rgrid.shape[0]*zgrid.shape[0])
    flphi = [0]*(rgrid.shape[0]*zgrid.shape[0])
    fli = 0

    # Interpolate field line shape from 100 to 1000 points, array values 
    # don't matter
    mxstep = 100
    coarse = np.linspace(0, 1, mxstep)
    fine = np.linspace(0, 1, 1000)
    for z in zgrid:
        for r in rgrid:
            # Trace mxstep steps, each 1e-2m with minus sign to go in 
            # counterclockwise direction
            fl = tracer.trace(r, z, mxstep=mxstep, ds=-1e-2, phistart=(316.5+40)*np.pi/180.)
            # Interpolation takes only first mxstep values; tracer sometimes 
            # gives more if field line hits wall
            flr[fli] = CubicSpline(coarse[:len(fl.R)], fl.R[:mxstep])(fine)
            flz[fli] = CubicSpline(coarse[:len(fl.R)], fl.Z[:mxstep])(fine)
            flphi[fli] = CubicSpline(coarse[:len(fl.R)], fl.phi[:mxstep])(fine)
            fli += 1

    return [flr, flz, flphi, rgrid, zgrid] 


if __name__ == '__main__':
    make_fls(1150611004, 0.78)
