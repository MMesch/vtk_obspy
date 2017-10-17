# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Filename: source.py
#  Purpose: Plots radiation patterns
# ---------------------------------------------------------------------

"""
Functions to compute and plot radiation patterns

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import numpy as np
from matplotlib.cm import get_cmap

from obspy.core.util import MATPLOTLIB_VERSION
from obspy.core.event.source import farfield
from obspy.imaging.scripts.mopad import MomentTensor, BeachBall
from obspy.imaging.mopad_wrapper import beach


def plot_radiation_pattern(mt, kind='mayavi', coordinate_system='RTP'):
    """
    Plot the P/S farfield radiation pattern on a unit sphere grid.

    The calculations are based on [Aki1980]_ eq. 4.29.

    :param mt: Focal mechanism NM x 6 (M11, M22, M33, M12, M13, M23 - the
        six independent components of the moment tensor, where the coordinate
        system is 1,2,3 = Up,South,East which equals r,theta,phi -
        Harvard/Global CMT convention). The relation to [Aki1980]_
        x,y,z equals North,East,Down convention is as follows: Mrr=Mzz,
        Mtt=Mxx, Mpp=Myy, Mrt=Mxz, Mrp=-Myz, Mtp=-Mxy.
    :param kind: One of:

        * **(A)** ``"mayavi"``: uses the mayavi library.
        * **(B)** ``"vtk"``: This vtk option writes two vtk files to the
          current working directory. ``rpattern.vtk`` contains the p and s
          wave farfield vector field. ``beachlines.vtk`` contains the nodal
          lines of the radiation pattern. A vtk glyph filter should be applied
          to the vector field (e.g. in ParaView) to visualize it.

    :type fig: :class:`matplotlib.figure.Figure`
    :param fig: Figure instance to use.
    :type show: bool
    :param show: Whether to show the figure after plotting or not. Can be
        used to do further customization of the plot before showing it.
    :returns: Matplotlib figure or ``None`` (if ``kind`` is ``"mayavi"`` or
        ``"vtk"``)
    """
    # reoorder all moment tensors to NED and RTP convention
    # name : COMPONENT              : NED sign and index
    # NED  : NN, EE, DD, NE, ND, ED : [0, 1, 2, 3, 4, 5]
    # USE  : UU, SS, EE, US, UE, SE : [1, 2, 0, -5, 3, -4]
    # RTP  : RR, TT, PP, RT, RP, TP : [1, 2, 0, -5, 3, -4]
    # DSE  : DD, SS, EE, DS, DE, SE : [1, 2, 0, -5, -3, 4]
    if coordinate_system == 'RTP' or coordinate_system == 'USE':
        signs = [1, 1, 1, -1, 1, -1]
        indices = [1, 2, 0, 5, 3, 4]
        ned_mt = [sign * mt[ind] for sign, ind in zip(signs, indices)]
    # the moment tensor has to be converted to
    # RTP/USE coordinates as well because the beachball routine relies
    # on it.
    # elif coordinate_system == 'DSE':
    #     signs = [1, 1, 1, -1, -1, 1]
    #     indices = [1, 2, 0, 5, 3, 4]
    #     ned_mt = [sign * mt[ind] for sign, ind in zip(signs, indices)]
    # elif coordinate_system == 'NED':
    #     ned_mt = mt
    else:
        msg = 'moment tensor in {:s} coordinates not implemented yet'
        raise NotImplementedError(msg.format(coordinate_system))

    if kind == 'mayavi':
        _plot_radiation_pattern_mayavi(ned_mt)

    elif kind == 'vtkfiles':
        # this saves two files, one with the vector field and one
        # with the nodal lines of the beachball
        fname_rpattern = 'rpattern.vtk'
        fname_beachlines = 'beachlines.vtk'
        _write_radiation_pattern_vtk(
            ned_mt, fname_rpattern=fname_rpattern,
            fname_beachlines=fname_beachlines)

    else:
        raise NotImplementedError('{:s} not implemented yet'.format(kind))


def _plot_radiation_pattern_mayavi(ned_mt):
    """
    Plot the radiation pattern using MayaVi.

    This private function uses the mayavi (vtk) library to plot the radiation
    pattern to screen. Note that you might have to set the QT_API environmental
    variable to e.g. export QT_API=pyqt that mayavi works properly.

    :param ned_mt: moment tensor in NED convention
    """
    # use mayavi if possible.
    try:
        from mayavi import mlab
    except Exception as err:
        print(err)
        msg = ("ObsPy failed to import MayaVi. "
               "You need to install the mayavi module "
               "(e.g. 'conda install mayavi', 'pip install mayavi'). "
               "If it is installed and still doesn't work, "
               "try setting the environmental variable QT_API to "
               "pyqt (e.g. export QT_API=pyqt) before running the "
               "code. Another option is to avoid mayavi and "
               "directly use kind='vtk' for vtk file output of the "
               "radiation pattern that can be used by external "
               "software like ParaView")
        raise ImportError(msg)

    # get mopad moment tensor
    mopad_mt = MomentTensor(ned_mt, system='NED')
    bb = BeachBall(mopad_mt, npoints=200)
    bb._setup_BB(unit_circle=False)

    # extract the coordinates of the nodal lines
    neg_nodalline = bb._nodalline_negative
    pos_nodalline = bb._nodalline_positive

    # add the first point to the end to close the nodal line
    neg_nodalline = np.hstack((neg_nodalline, neg_nodalline[:, 0][:, None]))
    pos_nodalline = np.hstack((pos_nodalline, pos_nodalline[:, 0][:, None]))

    # plot radiation pattern and nodal lines
    points = _equalarea_spherical_grid(nlat=20)
    dispp = farfield(ned_mt, points, type="P")
    disps = farfield(ned_mt, points, type="S")

    # get vector lengths
    normp = np.sum(dispp * points, axis=0)
    normp /= np.max(np.abs(normp))

    norms = np.sqrt(np.sum(disps * disps, axis=0))
    norms /= np.max(np.abs(norms))

    # make sphere to block view to the other side of the beachball
    rad = 0.8
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

    x = rad * sin(phi) * cos(theta)
    y = rad * sin(phi) * sin(theta)
    z = rad * cos(phi)

    # p wave radiation pattern
    mlab.figure(size=(800, 800), bgcolor=(0, 0, 0))
    pts1 = mlab.quiver3d(points[0], points[1], points[2],
                         dispp[0], dispp[1], dispp[2],
                         scalars=normp, vmin=-1., vmax=1.)
    pts1.glyph.color_mode = 'color_by_scalar'
    mlab.plot3d(*neg_nodalline, color=(0, 0.5, 0), tube_radius=0.01)
    mlab.plot3d(*pos_nodalline, color=(0, 0.5, 0), tube_radius=0.01)
    mlab.mesh(x, y, z, color=(0, 0, 0))

    # s wave radiation pattern
    mlab.figure(size=(800, 800), bgcolor=(0, 0, 0))
    pts2 = mlab.quiver3d(points[0], points[1], points[2],
                         disps[0], disps[1], disps[2], scalars=norms,
                         vmin=-0., vmax=1.)
    pts2.glyph.color_mode = 'color_by_scalar'
    mlab.plot3d(*neg_nodalline, color=(0, 0.5, 0), tube_radius=0.01)
    mlab.plot3d(*pos_nodalline, color=(0, 0.5, 0), tube_radius=0.01)
    mlab.mesh(x, y, z, color=(0, 0, 0))

    mlab.show()


def _write_radiation_pattern_vtk(
        ned_mt, fname_rpattern='rpattern.vtk',
        fname_beachlines='beachlines.vtk'):
    # output a vtkfile that can for exampled be displayed by ParaView
    mtensor = MomentTensor(ned_mt, system='NED')
    bb = BeachBall(mtensor, npoints=200)
    bb._setup_BB(unit_circle=False)

    # extract the coordinates of the nodal lines
    neg_nodalline = bb._nodalline_negative
    pos_nodalline = bb._nodalline_positive

    # plot radiation pattern and nodal lines
    points = _equalarea_spherical_grid()
    ndim, npoints = points.shape
    dispp = farfield(ned_mt, points, type="P")
    disps = farfield(ned_mt, points, type="S")

    # write vector field
    with open(fname_rpattern, 'w') as vtk_file:
        vtk_header = '# vtk DataFile Version 2.0\n' + \
                     'radiation pattern vector field\n' + \
                     'ASCII\n' + \
                     'DATASET UNSTRUCTURED_GRID\n' + \
                     'POINTS {:d} float\n'.format(npoints)

        vtk_file.write(vtk_header)
        # write point locations
        for x, y, z in np.transpose(points):
            vtk_file.write('{:.3e} {:.3e} {:.3e}\n'.format(x, y, z))
        # write vector field
        vtk_file.write('POINT_DATA {:d}\n'.format(npoints))
        vtk_file.write('VECTORS s_radiation float\n')
        for x, y, z in np.transpose(disps):
            vtk_file.write('{:.3e} {:.3e} {:.3e}\n'.format(x, y, z))
        vtk_file.write('VECTORS p_radiation float\n'.format(npoints))
        for x, y, z in np.transpose(dispp):
            vtk_file.write('{:.3e} {:.3e} {:.3e}\n'.format(x, y, z))

    # write nodal lines
    with open(fname_beachlines, 'w') as vtk_file:
        npts_neg = neg_nodalline.shape[1]
        npts_pos = pos_nodalline.shape[1]
        npts_tot = npts_neg + npts_pos
        vtk_header = '# vtk DataFile Version 2.0\n' + \
                     'beachball nodal lines\n' + \
                     'ASCII\n' + \
                     'DATASET UNSTRUCTURED_GRID\n' + \
                     'POINTS {:d} float\n'.format(npts_tot)

        vtk_file.write(vtk_header)
        # write point locations
        for x, y, z in np.transpose(neg_nodalline):
            vtk_file.write('{:.3e} {:.3e} {:.3e}\n'.format(x, y, z))
        for x, y, z in np.transpose(pos_nodalline):
            vtk_file.write('{:.3e} {:.3e} {:.3e}\n'.format(x, y, z))

        # write line segments
        vtk_file.write('\nCELLS 2 {:d}\n'.format(npts_tot + 4))

        ipoints = list(range(0, npts_neg)) + [0]
        vtk_file.write('{:d} '.format(npts_neg + 1))
        for ipoint in ipoints:
            if ipoint % 30 == 29:
                vtk_file.write('\n')
            vtk_file.write('{:d} '.format(ipoint))
        vtk_file.write('\n')

        ipoints = list(range(0, npts_pos)) + [0]
        vtk_file.write('{:d} '.format(npts_pos + 1))
        for ipoint in ipoints:
            if ipoint % 30 == 29:
                vtk_file.write('\n')
            vtk_file.write('{:d} '.format(ipoint + npts_neg))
        vtk_file.write('\n')

        # cell types. 4 means cell type is a poly_line
        vtk_file.write('\nCELL_TYPES 2\n')
        vtk_file.write('4\n4')


# ===== SUPPORT FUNCTIONS FOR SPHERICAL MESHES ETC STARTING HERE:
def _oriented_uv_sphere(ntheta=100, nphi=100, orientation=[0., 0., 1.]):
    """
    Returns a uv sphere (equidistant lat/lon grid) with its north-pole rotated
    to the input axis. It returns the spherical grid points that can be used to
    generate a QuadMesh on the sphere for surface plotting.

    :param nlat: number of latitudinal grid points (default = 100)
    :param nphi: number of longitudinal grid points (default = 100)
    :param orientation: axis of the north-pole of the sphere
                        (default = [0, 0, 1])
    """
    # make rotation matrix (after numpy mailing list)
    zaxis = np.array([0., 0., 1.])
    raxis = np.cross(orientation, zaxis)  # rotate z axis to null
    raxis_norm = np.linalg.norm(raxis)
    if raxis_norm < 1e-10:  # check for zero or 180 degree rotation
        rotmtx = np.eye(3, dtype=np.float64)
    else:
        raxis /= raxis_norm

        # angle between z and null
        angle = np.arccos(np.dot(zaxis, orientation))

        eye = np.eye(3, dtype=np.float64)
        raxis2 = np.outer(raxis, raxis)
        skew = np.array([[0, raxis[2], -raxis[1]],
                         [-raxis[2], 0, raxis[0]],
                         [raxis[1], -raxis[0], 0]])

        rotmtx = (raxis2 + np.cos(angle) * (eye - raxis2) +
                  np.sin(angle) * skew)

    # make uv sphere that is aligned with z-axis
    ntheta, nphi = 100, 100
    u = np.linspace(0, 2 * np.pi, nphi)
    v = np.linspace(0, np.pi, ntheta)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    # ravel point array and rotate them to the null axis
    points = np.vstack((x.flatten(), y.flatten(), z.flatten()))
    points = np.dot(rotmtx, points)
    return points


def _equalarea_spherical_grid(nlat=30):
    """
    Generates a simple spherical equalarea grid that adjust the number of
    longitude samples to the latitude. This grid is useful to plot vectors on
    the sphere but not surfaces.

    :param nlat: number of nodes in lat direction. The number of
                 nodes in lon direction is 2*nlat+1 at the equator
    """

    ndim = 3
    colats = np.linspace(0., np.pi, nlat)
    norms = np.sin(colats)
    # Scale number of point with latitude.
    nlons = (2 * nlat * norms + 1).astype(np.int_)

    # make colat/lon grid
    colatgrid, longrid = [], []
    for ilat in range(nlat):
        nlon = nlons[ilat]
        dlon = 2. * np.pi / nlon
        lons = np.linspace(0. + dlon / 2., 2. * np.pi - dlon / 2., nlon)
        for ilon in range(nlon):
            colatgrid.append(colats[ilat])
            longrid.append(lons[ilon])
    npoints = len(longrid)

    # get cartesian coordinates of spherical grid
    points = np.empty((ndim, npoints))
    points[0] = np.sin(colatgrid) * np.cos(longrid)
    points[1] = np.sin(colatgrid) * np.sin(longrid)
    points[2] = np.cos(colatgrid)

    return points
