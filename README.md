VTKOBS
======

Python toolbox for VTK visualizations of seismology related themes. This
package is based on obspy.

## Features
* 3D ray path visualizations
* vtk file export of ray paths
* moment tensor visualizations

## Example
this command:

```
./plot_rays.py --inv data/IU_stations.txt --phases P,Pdiff,PKP --evlat 0 --evlon 20
```

produces this plot:

![image](images/example1.png)


## Installation
The module should be downloaded, for example with:

```
git clone https://github.com/MMesch/vtkobs.git
```

and then executed from the installation directory
