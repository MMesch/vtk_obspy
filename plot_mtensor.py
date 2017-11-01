#!/usr/bin/env python

"""
Mayavi moment tensor visualizer.
"""

import sys
import os
import argparse

file_dir = os.path.dirname(os.path.realpath(__file__))
lib_dir = os.path.join(file_dir, 'lib')
sys.path.insert(0, lib_dir)

from vtk_obspy.source import plot_radiation_pattern  # NOQA

del sys.path[0]


def main():
    args = get_arguments()

    moment_tensor = [float(value) for value in args.mt.split(',')]
    if args.vtkfiles:
        kind = 'vtkfiles'
    else:
        kind = 'mayavi'
    # now plot rays
    plot_radiation_pattern(moment_tensor, kind=kind)


def get_arguments():
    """Return AttribDict with command line arguments."""
    description = ('3D visualization of moment tensor radiation '
                   'patterns.')
    usage = ('./plot_mtensor.py --mt 1,0,0,1,-1,0')

    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument('--mt', help='M11,M22,M33,M12,M13,M23')
    parser.add_argument('--vtkfiles', help='vtk files instead of mayavi',
                        action='store_true')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
