## Welcome to OpenMOOR

OpenMOOR is an open source cross-platform simulation program for static and dynamic analysis for mooring systems in offshore wind turbines and wave energy devices. 

## How to use
Compiled dynamic linking library files have been provided for use in Matlab for coupled analyis. Example file is also provided for use demonstration.
- MacOS
- Windows

## How to compile
We use [CMake](https://cmake.org) to setup and a MakeLists.txt is provided in the main folder along with the source files.
### MacOX
### Windows

## Dependencies
The following open source libraries are used by OpenMOOR:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix manipulation
- [rapidxml](http://rapidxml.sourceforge.net) for handling input file
- [boost/odeint](http://headmyshoulder.github.io/odeint-v2/) for integration

## View results in [Paraview](https://www.paraview.org)
### Prepare VTK files using Python

## Examples
[Cable model](https://github.com/chen-lin/OpenMOOR/blob/master/example/validation/Case3-5.gif?raw=true)
