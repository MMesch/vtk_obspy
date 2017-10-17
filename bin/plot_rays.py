#!/usr/bin/env python


"""
Mayavi ray path visualizer.
"""


import argparse
import obspy
from obspy import read_inventory, read_events

from ray_paths import plot_rays


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
        raise Exception('Specify station coordinates.')

    # now plot rays
    plot_rays(inventory, catalog)


def get_arguments():
    """Return AttribDict with command line arguments."""
    description = ('Lunch a 3D visualization of ray paths. The station and'
                   'event coordinates have to be given by the appropriate'
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

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
