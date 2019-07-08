Openmoor is an open source cross-platform simulation program for numerical simulation of statics and dynamics of mooring systems of offshore floating wind turbines and wave energy devices. In particular, it can consider the current of arbitrary profile.

### How to use
Openmoor has been compiled as dynamic linking libraries (DLL) for use in Matlab  environment. Check them in the _lib_ folder:

```
.
+-- lab
|   +-- MoorApiwin32.dll (DLL for win32 system, tested on windows 7)
|   +-- MoorApiwin64.dll (DLL for win64 system, tested on windows 7)
|   +-- libMoorApi.dylib (DLL for MacOS system, tested on Sierra)
|   +-- moorapi.h (Head file)
|   +-- CaseOC3.xml (Input file describing mooring system)
|   +-- current.dat (Input file describing current profile)
|   +-- openmoor_driver.m (Example Matlab code using openmoor DLL)
```
Set _lab_ as your work folder and then run openmoor_driver.m. For details on using openmoor and compiling openmoor on your own computer, please check the tutorial in pdf format in the _manual_ folder. 

The most updated information can be found on openmoor official [website](https://openmoor.github.io).

### A validation example
The scaled cable model described in this [paper](http://www.mdpi.com/2077-1312/4/1/5) is simulated. The unstretched length of the cable is 33 m. The computed upper end tensions agree well with the experimental data provided by the paper for these two cases.
- Case 1: forced upper end circular motion with radius 0.2 m and period 3.5 s
- Case 2: forced upper end circular motion with radius 0.2 m and period 1.25 s

Check reference [2] for the comparison of simulated responses and experimentally measured responses at the cable top end. The animations of the simulated cable motions in these two cases are shown below.

![](https://github.com/chen-lin/OpenMOOR/blob/master/examples/validation/Case3-5.gif?raw=true) | ![](https://github.com/chen-lin/OpenMOOR/blob/master/examples/validation/Case1-25.gif?raw=true)
:---:|:---:
Case 1 | Case 2

### Dependencies
The following open source libraries or third party functions are used by openmoor:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix manipulation
- [rapidxml](http://rapidxml.sourceforge.net) for handling input file
- [boost/odeint](http://headmyshoulder.github.io/odeint-v2/) for integration
- polyfit.cc

Check them and corresponding license in the _include_folder.

### Citation
Please cite [1,2] if you use openmoor in your own work. Other works relevant to openmoor is also listed below for your interest.

### References

[1] __Chen, L.__, Basu, B. & Nielsen, S.R.K. (2018). [A coupled finite difference mooring dynamics model for floating offshore wind turbine analysis](https://www.sciencedirect.com/science/article/pii/S0029801818307005). _Ocean Engineering,162_, 304-315.

[2] __Chen, L.__ & Basu, B. (2018). [Development of an open-source simulation tool for mooring systems](https://www.researchgate.net/publication/327424791_Development_of_an_open-source_simulation_tool_for_mooring_systems). In _Proceedings of the 2018 Civil Engineering Research in Ireland conference_ (CERI2018), Dublin, Ireland, pp. 823-828.

[3] __Chen, L.__ & Basu, B. (2019). [Wave-current interaction effects on structural responses of floating offshore wind turbines](https://onlinelibrary.wiley.com/doi/full/10.1002/we.2288). _Wind Energy, 22_(2), 327-339.

[4] __Chen, L.__, Basu, B. & Nielsen, S.R.K. (2019). [Nonlinear periodic response analysis of mooring cables using harmonic balance method](https://www.sciencedirect.com/science/article/pii/S0022460X18306126). _Journal of Sound and Vibration, 438_, 402-418.

[5] __Chen, L.__ & Basu, B. (2018). [Fatigue load estimation of a spar-type floating offshore wind turbine considering wave-current interactions](https://doi.org/10.1016/j.ijfatigue.2018.06.002). _International Journal of Fatigue, 116_, 421-428.

### License
Openmoor is licensed under the Apache License 2.0.

### Authors
- Dr. [Lin Chen](http://chen-lin.github.io)
- Prof. [Biswajit Basu](https://www.tcd.ie/research/profiles/?profile=basub)
