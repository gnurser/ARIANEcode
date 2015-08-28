from __future__ import print_function
from collections import OrderedDict
import fnmatch
import os, sys

from os.path import join as pjoin
from argparse import ArgumentParser


def extract(words):
    link = words[8]
    i = link.find('.nc')
    if i < 0:
        sys.exit('cannot find .nc in %s', line)
    link = link[i-5:i+3]
    filepath = words[-1]
    date_tail = filepath[-15:]
    return link, date_tail

if __name__ == '__main__':
    parser = ArgumentParser(description=
    """
    link fileset
    """)
    parser.add_argument('--linkfile','-f', dest='linkfile', default=None)
    parser.add_argument('--rootdir','-r', dest='rootdir', default=None)
    parser.add_argument('--linkdir', dest='linkdir', default=None)
    args = parser.parse_args()

    linkdict = OrderedDict()
    with open(args.linkfile) as f:
        lines = f.readlines()
        for line in lines:
            words = line.split()
            if len(words) > 9:
                link, date_tail = extract(words)
                linkdict[link] = date_tail

    for link, date_tail in linkdict.items():
	print(link, date_tail)

    matches = []
    for root, dirnames, filenames in os.walk(args.rootdir):
        for filename in fnmatch.filter(filenames, '*d05?.nc'):
            pathname = pjoin(root, filename)
            # print('found ', pathname)
            matches.append(pathname)

    for link, date_tail in linkdict.items():
        for match in matches:
            if date_tail in match:
                print('linking %s to %s' % (match, link))
                os.symlink(match, pjoin(args.linkdir,link))
