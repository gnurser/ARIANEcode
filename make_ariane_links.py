#!/usr/bin/env python
from __future__ import print_function
import sys, os
import numpy as np
from os.path import join as pjoin
import os.path
#from pathlib import Path
from argparse import ArgumentParser

def get_5d_files(rootdir,vble,
                  yrange=None, year0=1958, year1=2012):
    if yrange is None:
        yrange = range(year0, year1+1)
    tail = 'd05%s.nc' % vble
    print ('doing', vble)
    matches = []
    for year in yrange:
        year_str = '%4i' % year
        year_dir = pjoin(rootdir, year_str)
        print('%s' % year_str, end = ' ')
        year_files = os.listdir(year_dir)
        year_files.sort()
        matches += [pjoin(year_dir, file) for file in year_files if file.endswith(tail)]
    print()
    return matches

def adjust_links(vble, paths, skip_dates, dup_dates):
    if not (skip_dates or dup_dates):
        return
    for path in paths[:]:
        for skip in skip_dates[:]:
            if str(skip) in path:
                paths.remove(path)
                skip_dates.remove(skip)
                print('removed path ', path, 'included date', skip)
        for dup in dup_dates[:]:
            if str(dup) in path:
                insert_index = paths.index(path)
                paths.insert(insert_index, path)
                dup_dates.remove(dup)
                print('duplicated path ', path, 'duplicate date is', dup, ' inserted at ', insert_index)
        if not (skip_dates or dup_dates):
            return
    else:
        sys.exit('skip_dates', ' '.join(skip_dates),' and dup dates',' '.join(dup_dates),' not found')


def link_5d_files(rootdir, vbles = ('T','U','V','W'), link_dir = 'links',
                  yrange=None, year0=1958, year1=2012, test=False,
                   first_link=1, skip_dates=[], dup_dates=[]):
    if yrange is None:
        yrange = range(year0, year1+1)
    if not os.path.isdir(link_dir):
        os.makedirs(link_dir)
    for vble in vbles:
        paths = get_5d_files(rootdir, vble, yrange)
        print(skip_dates, dup_dates)
        adjust_links(vble, paths, skip_dates[:], dup_dates[:])
        for i,path in enumerate(paths):
            link = pjoin(link_dir,'%s%04i.nc' % (vble, i+first_link))
            if test:
                print(link, '  ',path)
            else:
                os.symlink(path, link)
        print(' *.%s-files :last link has index' % vble, len(paths) + first_link - 1)

if __name__ == '__main__':
    usage = "usage: %(prog)s oldstring newstring"
    parser = ArgumentParser(usage=usage)
    parser.add_argument('rootdir')
    parser.add_argument('-v','--variables',dest='variables',help='file suffices to link',nargs= '*', default=None)
    parser.add_argument('-y','--year_range',dest='year_range',help='first and last years to link', nargs= '*', type=int, default=None)
    parser.add_argument('-l','--link_dir',dest='link_dir',help='name of directory where links will go',default='links')
    parser.add_argument('-t','--test',dest='test', help="just test; don't make links", action='store_true', default=False)
    parser.add_argument('-f','--first',dest='first_link', help="number for first link",  type=int, default=1)
    parser.add_argument('-s','--skip',dest='skip_dates', help="dates to skip", nargs= '*',  type=int, default=[])
    parser.add_argument('-d','--dup',dest='dup_dates', help="dates to duplicate", nargs= '*',  type=int, default=[])
    args = parser.parse_args()

    if args.variables is None:
        args.variables = ('T','U','V','W')
    if args.year_range is None:
        args.year_range = (1958, 2012)
    year0, year1  = args.year_range
    link_5d_files(args.rootdir, vbles = args.variables, link_dir = args.link_dir,
                  yrange=None, year0=year0, year1=year1, test=args.test,
                   first_link=args.first_link, skip_dates=args.skip_dates, dup_dates = args.dup_dates)

