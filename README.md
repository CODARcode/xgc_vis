# README #

## Introduction ##

XGC_VIS provides tools to extract, track, and visualize topological blobs in XGC simulations.  

A set of tools are provided: 

* "Core": the core data analysis tool that can coule with the simulation via ADIOS in situ
* "GUI": the standalone Qt/OpenGL based GUI for visualizing analysis results in post-hoc

## Prerequisites ##

* [CMake](http://www.cmake.org/) (version >= 3.1.3)
* ADIOS (tested with version 1.11.0)
* MPI (tested with mpich-3.2)
* Boost (tested with )


## Build Guidelines ##

Please first clone the code from GitHub: 

`$ git clone git@github.com:CODARcode/xgc_vis.git`

Create a build directory for cmake: 

`$ mkdir build && cd build`


To build the code, replace "$your_adios_installation" with your ADIOS installation path in the following: 

``` shell
$ CC=mpicc CXX=mpicxx ADIOS_DIR=$your_adios_installation cmake ..
$ make
```

## Usage ##

```shell
$ ./core/xgcvis_core
./core/xgcvis_core:
  -m [ --mesh ] arg                  mesh_file
  -i [ --input ] arg                 input_file
  -o [ --output ] arg                output_file
  -r [ --read_method ] arg (=BP)     read_method (BP|DATASPACES|DIMES|FLEXPATH)
  -w [ --write_method ] arg (=POSIX) write_method (POSIX|MPI)
  --h                                display this information
```

For example

```shell
$ ./core/xgcvis_core ~/workspace/data/xgc_feb2018/xgc.mesh.bp ~/workspace/data/xgc_feb2018/xgc.3d.00450.bp output.bp --read_method=BP --write_method=POSIX
==========================================
filename_mesh=/Users/hguo/workspace/data/xgc_feb2018/xgc.mesh.bp
filename_input=/Users/hguo/workspace/data/xgc_feb2018/xgc.3d.00450.bp
filename_output=output.bp
read_method=BP
write_method=POSIX
==========================================
reading mesh...
nNodes=125098, nTriangles=248866
reading data...
starting analysis..
building contour tree for plane 0, nNodes=125098
dumping results..
done.
```

Then you can open output.bp with VisIt.


## Contact ##

* [Hanqi Guo](http://www.mcs.anl.gov/~hguo/), [hguo@anl.gov](mailto:hguo@anl.gov)
* Tom Peterka, [tpeterka@mcs.anl.gov](mailto:tpeterka@mcs.anl.gov)
