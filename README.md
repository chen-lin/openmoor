## Welcome to OpenMOOR

OpenMOOR is an open source cross-platform simulation program for static and dynamic analysis for mooring systems in offshore wind turbines and wave energy devices. 

## How to use
OpenMOOR can be used as a standalone program or more frequently as dynamic linking library for coupled analysis. 
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
### Validation
The scaled model cable described in [paper](http://www.mdpi.com/2077-1312/4/1/5) is simulated.
- Case 1: upper end forced motion with period 3.5 s
![case1](https://github.com/chen-lin/OpenMOOR/blob/master/example/validation/Case3-5.gif)

- Case 2: upper end forced motion with period 1.25 s
<img src="https://github.com/chen-lin/OpenMOOR/blob/master/example/validation/Case1-25.gif" width="500">
