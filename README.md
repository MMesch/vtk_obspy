VTKOBS
======


Python toolbox for VTK visualizations of seismology related themes. This
package is based on obspy.


Features
--------

* 3D ray path visualizations
* vtk file export of ray paths


Installation
------------
The module just needs to be downloaded, for example with:
```
git clone https://github.com/MMesch/vtkobs.git
```
The plot scripts can be executed directly from the installation directory.
Another option is to install the python routines to the system with
```
pip install .
```
before using the plot scripts.

### interactive plotting
For interactive plotting you need to install `mayavi 4.5.0` and make sure that
it works. In anaconda this can be done with:
```
conda install -c menpo mayavi
```
Be careful, you might need to downgrade pyqt to version 4 instead of 5.


Example
-------

### VTK file export
plotting with the ``--vtkfiles`` option produces vtk files that can be
visualized with paraview.

```
./plot_rays.py --inv data/IU_stations.txt --phases P,Pdiff,PKP --evlat 0
               --evlon 20 --vtkfiles
```

### interactive plot
this command:

```
./plot_rays.py --inv data/IU_stations.txt --phases P,Pdiff,PKP --evlat 0 --evlon 20
```

produces this plot:

![image](images/example1.png)
