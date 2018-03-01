## Welcome to OpenMOOR

OpenMOOR is an open source cross-platform simulation program for numerical simulation of statics and dynamics of mooring systems of offshore floating wind turbines and wave energy devices. 

## Demonstration
### Validation
The scaled model cable of unstretched length of 33 m tested in this [paper](http://www.mdpi.com/2077-1312/4/1/5) is simulated. The computed upper end tensions agree well with the experimental data provided by the paper for these two cases.
- Case 1: forced upper end motion with radius 0.2 m and period 3.5 s
![case1](https://github.com/chen-lin/OpenMOOR/blob/master/example/validation/Case3-5.gif?raw=true)

- Case 2: forced upper end motion with radius 0.2 m and period 1.25 s
![case2](https://github.com/chen-lin/OpenMOOR/blob/master/example/validation/Case1-25.gif?raw=true)

### Coupled with FOWT

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


