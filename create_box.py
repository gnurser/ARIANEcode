#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import os, sys
from os.path import join as pjoin
from argparse import ArgumentParser, ArgumentTypeError
import datetime


import numpy as np
import numpy.ma as ma
import netCDF4
from netCDF4 import Dataset

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from traj_movie import get_circles

def create_initial_positions(xdivide=10, ydivide=10,lon0=55.,lon1=56.,lat0=-21.5, lat1=-20.5
                            ,runlength=147, kdepth=-1.5, outfile='initial_positions.txt', res='0083', crop=False):

    lons = np.linspace(lon0, lon1, xdivide+1)
    lats = np.linspace(lat0, lat1, ydivide+1)

    domaindir=pjoin('/Users/agn/Data/NEMO',res)
    path = pjoin(domaindir,'mask.nc')
    print('path for mask file', path)
    with Dataset(path) as f:
        Ndlat = f.variables['nav_lat']
        meridional = Ndlat[:,0]
        jeq = meridional.searchsorted(0.)
        Ndlon = f.variables['nav_lon']
        zonal = Ndlon[jeq,:]

        ibreak = np.argmax(zonal) + 1
        if lon0 > zonal[0] and lon1 < zonal[ibreak-1]:
            zonal_part = zonal[:ibreak]
            iadd = 0
        elif lon0 > zonal[ibreak] and lon1 < zonal[-1]:
            zonal_part = zonal[ibreak:]
            iadd = ibreal
        else:
            sys.exit('lon0 and lon1 bracket dateline; not implemented yet')

        i0 = zonal_part.searchsorted(lon0) - 1 + iadd
        i1 = zonal_part.searchsorted(lon1) + 1 + iadd
        imid = (i0 + i1)//2
        meridional = Ndlat[:,imid]
        j0 = meridional.searchsorted(lat0) - 1
        j1 = meridional.searchsorted(lat1) + 1
        jmid = (j0 + j1)//2
        zonal = Ndlon[jmid,:]

        dlon = np.diff(zonal)
        ibreak = np.argmin(dlon)

        ri = np.interp(lons, zonal[i0:i1], np.arange(i0,i1,dtype=zonal.dtype))
        rj = np.interp(lats, meridional[j0:j1], np.arange(j0,j1,dtype=meridional.dtype))
        print('specified longitudes are','\n',lons)
        print('i values are','\n',ri + 1. -.5)
        print('specified latitudes are','\n',lats)
        print('j values are','\n',rj + 1. -.5)

        i00, i11 = i0 - 1, i1 + 1
        j00, j11 = j0 - 1, j1 + 1
        tmask = ~(f.variables['tmaskutil'][0,j00:j11,i00:i11].astype(np.bool))

    uvmask = tmask[1:-1,1:-1] + tmask[:-2,1:-1] + tmask[2:,1:-1] \
                + tmask[1:-1,:-2] + tmask[1:-1,2:]

    if crop:
        inner_circle, outer_circle, fuel_circle, box = get_circles()
        polygon = Polygon( zip(*box) )
    npoints = 0
    with open(outfile, 'w') as f:
        for y, lat in zip(rj, lats):
            for x, lon in zip(ri, lons):
                if (not crop) or (crop and polygon.contains(Point(lon, lat))):
                    i, j = int(x), int(y)
                    if not uvmask[j-j0, i-i0]:
                        # ri and rj are  c-indices relative to T-cells
                        # subtract .5 to produce u v indices as required by ariane
                        # ... and add 1 to convert from C to fortran numbering
                        f.write('%10g%10g%10g%10g%5g\n' %(x + 1. -.5, y + 1. -.5 , kdepth, runlength, 1))
                        npoints +=1
    print('# of sea points started is %g out of %g' % (npoints, len(rj)*len(ri)))

if __name__ == '__main__':
    parser = ArgumentParser(description=
    """
    produce movie of particle positions
    """)
    parser.add_argument('--res',dest='res',
                        help="model resolution ",default='0083')
    parser.add_argument('--nx',dest='xdivide',
                        help="x divisions ",type=int, default=100)
    parser.add_argument('--ny',dest='ydivide',
                        help="y divisions ",type=int, default=100)
    parser.add_argument('--lon0',dest='lon0',
                        help="W bdry ",type=float, default=55.)
    parser.add_argument('--lon1',dest='lon1',
                        help="E bdry ",type=float, default=56.)
    parser.add_argument('--lat0',dest='lat0',
                        help="S bdry ",type=float, default=-21.5)
    parser.add_argument('--lat1',dest='lat1',
                        help="N bdry ",type=float, default=-20.5)
    parser.add_argument('--time_index',dest='time_index',
                        help="time index when particles at initial positions ",type=int, default=1)
    parser.add_argument('--kdepth',dest='kdepth',
                        help="depth in terms of k (k=1 is ths surface) ",type=float, default=-1.5)
    parser.add_argument('--outfile', '-o',dest='outfile',
                        help="output file of points ",default='initial_positions.txt')
    parser.add_argument('--crop', '-c',dest='crop', action='store_true',
                        help="crop by satellite zone ",default=False)
    args = parser.parse_args()

    create_initial_positions(xdivide=args.xdivide, ydivide=args.ydivide,lon0=args.lon0,lon1=args.lon1,lat0=args.lat0, lat1=args.lat1
                            ,runlength=args.runlength, kdepth=args.kdepth, outfile=args.outfile, res=args.res, crop=args.crop)
