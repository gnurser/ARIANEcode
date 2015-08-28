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

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path as MplPath
from matplotlib.patches import Polygon

import small_circle


class Savefiles(object):
    def __init__(self,movie_name=None,pdf=False,dpi=None,tempdir='temp',fmt='%04i'):
        self.movie_name = movie_name
        self.fmt=fmt
        if movie_name is None and pdf:
            self.suffix = 'pdf'
        else:
            self.suffix = 'png'

        if dpi==None:
            if movie_name is not None:
                self.dpi=100
            else:
                self.dpi=300
        else:
            self.dpi = dpi

        if  movie_name is not None:
            tempdir += movie_name
            if os.path.exists(tempdir):
                oldfiles = os.listdir(tempdir)
                [os.remove(pjoin(tempdir,file)) for file in oldfiles]
            else:
               os.mkdir(tempdir)

            self.files = []
            self.tempdir = tempdir
        else:
            self.tempdir = '.'

    def savefile(self,fig,nstep=None,fname='_tmp'):
        if nstep is not None:
            fname += (self.fmt+'.') %nstep + self.suffix
        print('Saving file', fname)
        fig.savefig(pjoin(self.tempdir,fname),dpi=self.dpi)
        del fig
        self.files += [fname]
        return fname

    def make_movie(self,clean=True, use_ffmpeg=True):
        print('Making movie animation.mpg - this make take a while')
        currdir = os.getcwd()
        os.chdir(self.tempdir)
        movie_path = pjoin(currdir,self.movie_name)
        if use_ffmpeg:
            command = 'ffmpeg -i _tmp0%3d.png -filter:v "setpts=2.0*PTS" -f mp4 -vcodec h264 -pix_fmt yuv420p {0}.mp4'.format(movie_path)
        else:
            command = 'mencoder -nosound -ovc x264 -x264encopts bitrate=1200:threads=2 -mf type=png:fps=30 -o {0}.avi mf://\*.png -v'.format(movie_path)
        fail = os.system(command)
        # print(fail, flush=True)
        # cleanup
        if clean and fail==0:
            for fname in self.files: os.remove(fname)
            os.chdir(currdir)
            os.rmdir(self.tempdir)


class TRAJ(object):
    @classmethod
    def set_root(cls, rootdir='.'):
        """
        Set rootdir
        """
        TRAJ.root = rootdir

    def __init__(self, run,  color='b', rootdir=None):
        self.color = color
        if run[1] == 'p':
            self.name = '%s %% wind' % run[0]
        else:
            self.name = run
        if rootdir is None:
            rootdir = TRAJ.root
        rundir=pjoin(rootdir,run)
        pathname = pjoin(rundir,'ariane_trajectories_qualitative.nc')
        if not os.path.exists(pathname):
            sys.exit('cannot find path %s' % pathname)
        with Dataset(pathname) as f:
            self.sense = f.forback
            fv = f.variables
            self.init_x = fv['init_x'][:]
            self.init_y = fv['init_y'][:]
            self.init_z = fv['init_z'][:]
            self.init_s = fv['init_age'][:]
            self.traj_lon, self.traj_lat, self.traj_depth, self.traj_time = \
                [fv[x][...] for x in ('traj_lon', 'traj_lat', 'traj_depth', 'traj_time')]
            self.final_z = fv['final_z'][:]

    def subsample(self, vertices):
        xiyi = np.vstack((self.traj_lon[0,:], self.traj_lat[0,:])).T
        mpath = MplPath( vertices )
        outside_area = ~mpath.contains_points(xiyi)
        self.init_x, self.init_y, self.init_z, self.init_s, \
                   self.traj_lon, self.traj_lat, self.traj_depth, self.traj_time, self.final_z = \
                 [x.compress(outside_area,axis=-1)  for x in (self.init_x, self.init_y, self.init_z, self.init_s,
                   self.traj_lon, self.traj_lat, self.traj_depth, self.traj_time, self.final_z)]

    def ptraj(self, n):
        print(self.traj_lon[:,n], self.traj_lat[:,n], self.traj_depth[:,n], self.traj_time[:,n])
        print(1*'\n')


    def ptraj_start(self):
        print(self.traj_lon[0,:], 2*'\n', self.traj_lat[0,:], 2*'\n', self.traj_depth[0,:])


    def plot(self, n, ax=None, m=None, color=None):
        if color is None:
            color = self.color

        if m is not None:
            graph = m
            x, y = m(self.traj_lon[:,n], self.traj_lat[:,n])
        else:
            if ax is None:
                try:
                    ax = plt.gca()
                except:
                    fig, ax = plt.subplots()

            graph = ax
            x, y = self.traj_lon[:,n], self.traj_lat[:,n]

        graph.plot(x, y, c=color)

    def plot_dots(self, nt=100, ax=None, m=None, color=None, **kwargs):
        if color is None:
            color = self.color

        if m is not None:
            graph = m
            x, y = m(self.traj_lon[nt,:], self.traj_lat[nt,:])
        else:
            if ax is None:
                try:
                    ax = plt.gca()
                except:
                    fig, ax = plt.subplots()

            graph = ax
            x, y = self.traj_lon[nt,:], self.traj_lat[nt,:]

        return graph.scatter(x, y, lw=0, c=color, **kwargs)
        # realdots = graph.scatter(x, y, lw=0, c=color, **kwargs)
        # return graph.plot([np.NaN],[np.NaN], 'o', c=color, lw=0, alpha=0.7, markersize=15)


def south_indian():
    fig, ax = plt.subplots()
    # m = Basemap(projection='merc',llcrnrlat=-45,urcrnrlat=10,
    m = Basemap(projection='cea',llcrnrlat=-45,urcrnrlat=10,
                    llcrnrlon=30,urcrnrlon=130,lat_ts=20,resolution='l',ax=ax)
    parallels = np.linspace(-50.,10.,7)
    mlabels = [0,0,0,1]
    plabels = [1,0,0,0]
    meridians = np.linspace(-180.,160.,18)

    xlines = {'linewidth':0.3,'dashes':[1,0.1],'size':'x-small'}
    m.drawparallels(parallels,labels=plabels,**xlines)
    m.drawmeridians(meridians,labels=mlabels,**xlines)
    conts = m.fillcontinents(color='k',lake_color='w')

    return fig, ax, m

def plot_trajs(trajlist, runlist, savefig=True):
    fig, ax, m = south_indian()

    for ntraj in trajlist:
        if min([U.traj_time[:,ntraj].max() for U in runlist])> .99:
            for U in runlist:
              U.plot(ntraj, m=m)

    if savefig:
        fig.savefig('trajs.png', dpi=300)


class SearchArea(object):
    def __init__(self, lat0=-38.5, lat1=-33.5):
        self.lat0, self.lat1 = lat0, lat1
        self.lat00, self.lat11 = lat0 -.5, lat1 + .5

    def get_sw_indian(self, lon_lat):
        lon, lat = lon_lat
        ok = (lon > 80.) * (lat > self.lat00) * (lat < self.lat11)
        clon, clat =  [ma.masked_where(~ok, x).compressed() for x in (lon, lat)]
        # ensure clat is increasing
        if clat[1]<clat[0]:
            clat = clat[::-1]
            clon = clon[::-1]

        return clon, clat

    def box(self, xy1, xy2):
        # print(xy1)
        # print(xy2)
        ndivs = 20
        npts = ndivs + 1
        lats_inner = np.linspace(self.lat0, self.lat1, npts)
        lats = np.empty([2*npts])
        lats[:npts] = lats_inner
        lats[npts:] = lats_inner[::-1]
        # print(lats)
        lons = np.empty_like(lats)
        lons_inner = np.interp(lats[:npts], xy1[1], xy1[0])
        lons_outer = np.interp(lats[npts:], xy2[1], xy2[0])
        lons[:npts] = (3.*lons_inner + 1.*lons_outer[::-1])/4.
        lons[npts:] = (3.*lons_outer + 1.*lons_inner[::-1])/4.
        # print(lons)
        return lons, lats

def get_circles():
    radius = 6731.*1.e3
    meters = 1.e4

    # as given in http://www.earth.com/2014/03/flight-mh370-search-data-in-google-earth/
    # igrex_center = (94. + 25./60., 9.75)
    # correcting longitude so consistent with ECMWF
    igrex_center = (97.1, 9.75)

    # taken by eye from ECMWF plot... southernmost extent
    igrex_point = igrex_center[0], -40.25

    # position of Inmarsat-3 F1
    sat_center = 64.5, 0.
    # estimated from intersections of ECMWF lines with Javanese coastline, using http://www.latlong.net
    satpoint_inner = 105.3, -6.78
    satpoint_outer = 108.56, -6.75

    clon, clat = sat_center
    satellite_circles = []
    for point in  (satpoint_inner, satpoint_outer):
        llon, llat = point
        satellite = small_circle.small_circle(radius,clon,clat,llon,llat,meters)
        satellite_circles.append(satellite)

    clon, clat = igrex_center
    llon, llat = igrex_point
    fuel_circle = small_circle.small_circle(radius,clon,clat,llon,llat,meters)

    # get arrays each of shape (2,npoints) ...
    inner_circle, outer_circle = [np.array(x).T.copy() for x in satellite_circles]
    fuel_circle = np.array(fuel_circle).T.copy()

    # area of maximum likelihood
    # from http://jacc.gov.au/media/reports/2014/october/files/
    #  AE-2014-054_MH370-FlightPathAnalysisUpdate.pdf
    area = SearchArea(lat0=-38.5, lat1=-33.5)
    inner_sw_indian, outer_sw_indian = [area.get_sw_indian(xy) for xy in (inner_circle, outer_circle)]
    box = area.box(inner_sw_indian, outer_sw_indian)

    return inner_circle, outer_circle, fuel_circle, box


def plot_circles(m=None):
    inner_circle, outer_circle, fuel_circle, box = get_circles()
    for circle in inner_circle, outer_circle:
        lons, lats = circle
        x, y = m(lons,lats)
        m.plot(x, y, color='k', lw=0.7)

    lons, lats = fuel_circle
    x, y = m(lons,lats)
    m.plot(x, y, color='gray', lw=1.)

    lons, lats = box
    x, y = m(lons,lats)
    poly = Polygon(zip(x, y), facecolor='m', alpha=0.5, edgecolor='none')
    m.ax.add_patch(poly)


def leapday_after(plot_date, end_date, leapyears=[2004, 2008, 2012]):
    for leapyear in leapyears:
        leap_date = datetime.datetime(leapyear, 2, 29)
        if plot_date<leap_date and end_date>leap_date:
            return 1
    else:
        return 0

def do_dots(l, runs, savefile=None, reverse=True, circles=True, ntmax=147,
             dotsize=20, alpha=1., date=True):
    fig, ax, m = south_indian()
    if reverse:
        nt=ntmax - l - 1
    else:
        nt = l
    print('l=',l,'nt=',nt)

    if circles:
        plot_circles(m=m)

    leglist = []
    for U in runs:
        U.plot_dots(nt=nt, m=m, s=dotsize, alpha=alpha)
        legdot = mlines.Line2D([], [], color=U.color, marker='.',
                          markersize=8, lw=0, label=U.name, alpha = 0.6)
        leglist.append(legdot)

    ax.legend(handles=leglist, loc='upper left')

    if date:
        if U.sense == 'backward':
            plot_date = U.end_date - datetime.timedelta(days = nt*5)
        else:
            plot_date = U.crash_date + datetime.timedelta(days = nt*5)

        plot_date -= datetime.timedelta(days=leapday_after(plot_date, U.end_date))
        print(U.end_date.strftime('%d %b %Y'))
        print(plot_date.strftime('%d %b %Y'))
        title_text = plot_date.strftime('%d %b %Y')
    else:
        if U.sense == 'backward':
            title_text = ' %g days of backtracking from Reunion' % (5*nt,)
        else:
            title_text = ' %g days of tracking from crash site' % (5*nt,)

    ax.set_title(title_text)

    if savefile is not None:
        savefile.savefile(fig, nstep=l)
    else:
        fig.savefig('dots%03i.png' % l, dpi=300)
        fig.savefig('dots%03i.pdf' % l)

    try:
        plt.close(fig)
    except:
        pass


def valid_date(s):
    """
    Checks that input dates are in the right format, and
    read them into datetime objects
    """
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise ArgumentTypeError(msg)

# class RUNS(object):
#     pass

if __name__ == '__main__':
    parser = ArgumentParser(description=
    """
    produce movie of particle positions
    """)
    parser.add_argument('--movie','-m',dest='movie',
                        help='name of movie',
                        nargs='?', const='MH370_dots', default=None)
    parser.add_argument('--keep_pngs','-k',dest='keep_pngs',
                        help="don't clean out pngs",action='store_true',
                        default=False)
    parser.add_argument('-r',dest='rootdir',
                        help="directory of ARIANE runs ",default='.')
    parser.add_argument('--ntmax',dest='ntmax',
                        help="max time length of trajectories ",type=int, default=None)
    parser.add_argument('--dotsize',dest='dotsize',
                        help="area of dots in points^2 ",type=int, default=20)
    parser.add_argument('--dotalpha',dest='dotalpha',
                        help="transparency of dots ",type=float, default=0.8)
    parser.add_argument('--dpi', dest='dpi',
                        help="resolution of pngs in dpi ",type=float, default=300.)
    parser.add_argument(dest='runs',
                        help="list of runs to plot ", nargs='+', default=None)
    parser.add_argument('--date', dest='date',
                        help="put date of 5 day period instead of day number on plots ",
                        action='store_true', default=False)
    parser.add_argument('--subarea', dest='subarea',
                        help="vertices of valid area to subsample initial"
                        " positions of particles, as successive lon-lat pairs",
                        type=float, nargs='+', default=None)
    parser.add_argument('--end_date',dest='end_date',
                        help="The End Date format YYYY-MM-DD ",type=valid_date,
                         default="2010-07-29")
    parser.add_argument('--crash_date',dest='crash_date',
                        help="The crash date format YYYY-MM-DD ",type=valid_date,
                         default="2009-03-08")
    args = parser.parse_args()

    clean_pngs = not args.keep_pngs


    TRAJ.set_root(args.rootdir)
    TRAJ.end_date = args.end_date
    TRAJ.crash_date = args.crash_date

    nsea_days = args.end_date - args.crash_date
    lcrash = nsea_days.days//5
    model_crash_date = args.end_date - datetime.timedelta(days=5*lcrash)
    model_crash_date -= datetime.timedelta(days=leapday_after(model_crash_date, args.end_date))
    print('model crash date is %s' % model_crash_date.strftime('%d %b %Y'))
    colors = ['b','r','g','c']
    ulist = []
    for i,run in enumerate(args.runs):
        ulist.append(TRAJ(run, color=colors[i]))

    ntmax, ntrajs = ulist[0].traj_lon.shape
    if args.ntmax is not None:
        ntmax = args.ntmax

    print('number of trajectories =', ntrajs)
    print('length of trajectories =', ntmax)
    
    if args.subarea is not None:
        vertices = np.array(args.subarea).reshape(-1,2)
        print('\nOnly keeping points in polygon bound by ',
              [tuple(vertices[i]) for i in range(vertices.shape[0])])
        for U in ulist:
            U.subsample(vertices)
        print('  number of remaining trajectories is ',ulist[0].init_x.shape[0])


    # processes = max(cpu_count() - 2,1)
    # sys.stdout.flush()
    # if processes>1:
    #     # create processes workers
    #     pool = Pool(processes=processes)
    #     xyzs_list = pool.map(dopict,self.traj_list)
    #     # wait for all workers to finish & exit parallellized stuff
    #     pool.close()
    # else:
    #     xyzs_list = []
    #     for traj in self.traj_list:
    #         xyzs_list.append( dopict(traj) )
    # print()


    do_dots(lcrash, ulist, savefile=None, reverse=False,
             date=args.date, dotsize=args.dotsize, alpha=args.dotalpha)

    if args.movie:
        save = Savefiles(movie_name=args.movie, dpi = args.dpi)
        for l in range(ntmax):
            do_dots(l, ulist, savefile=save, reverse=False,
                    date=args.date, dotsize=args.dotsize, alpha=args.dotalpha)
        save.make_movie(clean=clean_pngs)
