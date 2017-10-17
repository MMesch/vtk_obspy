#!/usr/bin/env python


"""
Mayavi ray path visualizer.
"""


import argparse
import obspy
from obspy import read_inventory, read_events

from vtk_obspy.ray_paths import plot_rays


def main():
    args = get_arguments()

    # fill catalog and inventory objects that are used as input for the
    # ray path plotting routine
    if not (args.evlat is None and args.evlon is None):
        otime = obspy.UTCDateTime()
        origin = obspy.core.event.Origin(latitude=float(args.evlat),
                                         longitude=float(args.evlon),
                                         depth=float(args.evdep),
                                         time=otime)
        magnitude = obspy.core.event.Magnitude(mag=7.)
        event = obspy.core.event.Event(origins=[origin],
                                       magnitudes=[magnitude])
        catalog = obspy.core.event.Catalog(events=[event])
    elif args.cat:
        catalog = read_events(args.cat)
    else:
        raise Exception('Please specify the event coordinates.')
    if not (args.stlat is None and args.stlon is None):
        station = obspy.core.inventory.Station(
                    code='STA', latitude=float(args.stlat),
                    longitude=float(args.stlon), elevation=0.)
        network = obspy.core.inventory.Network(
                    code='NET', stations=[station])
        inventory = obspy.core.inventory.Inventory(
                source='ME', networks=[network])
    elif args.inv:
        inventory = read_inventory(args.inv)
    else:
        raise Exception('Please specify the station coordinates.')

    phase_list = [phase for phase in args.phases.split(',')]
    if args.vtkfiles:
        kind = 'vtkfiles'
    else:
        kind = 'mayavi'
    # now plot rays
    plot_rays(inventory, catalog, phase_list=phase_list,
              station_labels=not args.nostalabels,
              event_labels=not args.noevlabels, kind=kind,
              colorscheme=args.colorscheme)


def get_arguments():
    """Return AttribDict with command line arguments."""
    description = ('3D visualization of ray paths. The station and '
                   'event coordinates have to be given by the appropriate '
                   'arguments.')
    usage = ('./plot_rays.py --cat catalog.xml --inv inventory.xml\n'
             './plot_rays.py --cat catalog.xml --stlat 20 --stlon 30')

    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument('--cat', help='path to catalog file')
    parser.add_argument('--inv', help='path to inventory file')
    parser.add_argument('--evlat', help='event latitude')
    parser.add_argument('--evlon', help='event longitude')
    parser.add_argument('--evdep', help='event depth', default=0.0)
    parser.add_argument('--stlat', help='station latitude')
    parser.add_argument('--stlon', help='station longitude')
    parser.add_argument('--phases', help='e.g. P,PP', default='P')
    parser.add_argument('--noevlabels', help='switch off event labels',
                        action='store_true')
    parser.add_argument('--nostalabels', help='switch off station labels',
                        action='store_true')
    parser.add_argument('--vtkfiles', help='vtk files instead of mayavi',
                        action='store_true')
    parser.add_argument('--colorscheme', help='dark or bright', default='dark')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
