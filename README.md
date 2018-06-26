[![Build Status](https://travis-ci.org/ALPSCore/CT-INT.svg?branch=master)](https://travis-ci.org/ALPSCore/CT-INT)

CT-INT
======

Interaction expansion quantum Monte Carlo impurity solver: An open-source implementation of the continuous-time interaction-expansion quantum Monte Carlo method.

This program solves impurity models with onsite Coulomb interactions and complex Weiss functions.
The code is built on the ALPSCore libraries (Applications and Libraries for Physics Simulations Core libraries).

The code was written by Hiroshi Shinaoka in collaboration with Yusuke Nomura based on its ancestor written by Emanuel Gull.
All the files expect for spline.h are licensed under GPLv3 or later.

# Table of Contents
- [Requirements](#requirements)
- [Manual source installation](#manual-source-installation)
- [Trouble shooting](#trouble-shooting)

## Requirements
### C++11 compiler
Modern C++ compier supporting C++11

### ALPSCore libraries (>= 2.1)
The ALPSCore libraries are need to be properly installed, see [ALPSCore library](https://github.com/ALPSCore/ALPSCore). 
If your installation of ALPSCore was built with C++03 or C++98, please rebuild it with C++11 or C++14. See [wiki](https://github.com/ALPSCore/ALPSCore/wiki/Installation#manual-source-installation)

### Boost (>= 1.54.0)
Only header-file libraries are needed. The dependencies will be taken care of by ALPSCore.

### Eigen3
The latest ALPSCore libraries depend on Eigen3.
Thus, this dependency will automatically propagate to the build of the CT-INT.

### MPI
MPI environment is required.

## Manual source installation
The solver depends on ALPSCore libraries and some Boost header-file only libraries.
These libraries must be preinstalled.
The CT-INT does NOT depend on any Boost binary libraries.

The CT-INT package can be obtained by following methods:
* Clone Git repository at Github
```
$ git clone https://github.com/ALPSCore/CT-INT.git
```

Then, make a (separated) build directly, and provide something like:

```
$ mkdir build
$ cd build
$ cmake\
$     -DALPSCore_DIR=/path/to/ALPSCore \
$     -DCMAKE_INSTALL_PREFIX=/path/to/install/dir \
$     -DCMAKE_CXX_COMPILER=/path/to/C++/compiler \
$    ../CT-INT
$ make
$ make test
$ make install
```
You must use a MPI C++ compiler.
This may be done by setting the path of your MPI wrapper compiler (such as mpic++) to CMAKE\_CXX\_COMPILER.
Note that, in such cases, MPI must be enabled also in the installation of ALPSCore.
If cmake does not find boost, please tell cmake the installation directory of boost by using the option "-DBOOST_ROOT=***".
If cmake does not find Eigen3, please set use the option "-DEIGEN3\_INCLUDE\_DIR".
For instance, if the Core file of your Eigen3 is located at "/opt/local/include/eigen3/Eigen/Core",
please use "-DEIGEN3\_INCLUDE\_DIR=/opt/local/include/eigen3".

Please make sure that ALPSCore/CT-INT is going to be built with the same C++ standard
as that used for building the ALPSCore libraries (>= C++11).


## Trouble shooting
* Some libraries are not found at runtime.<br>
When you install the executalbe to your installation path by "make install", CMake removes the paths of dynamic libraries from the binary.
When you launch "/path/to/install/dir/hybmat", some dynamic libraries which were visible in the build may not be found.
In this case, please set your environment variables correctly (e.g., LD\_LIBRARY\_PATH) so that the system can find these libraries at runtime. More information is found [here]
(https://cmake.org/Wiki/CMake_RPATH_handling).