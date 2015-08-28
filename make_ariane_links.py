import os
from os.path import join as pjoin
max_n = 2263
for suffix in ('U','V','W','T'):
    filepath = pjoin(os.getcwd(), 'ORCA0083-N01_20100709d05%s.nc' % suffix)
    for n in range(max_n):
        lname = '%s%04i.nc' % (suffix, n)
        symlink = pjoin('../test_means_5d', lname)
        os.symlink(filepath, symlink)
