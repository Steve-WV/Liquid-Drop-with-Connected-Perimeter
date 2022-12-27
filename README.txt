Guidelines how to use the code: 

REQUIREMENTS:
 - Cmake version no less than 3.7.2
 - A suitable C++ compiler
 - The linear algebra library Eigen
 - The Boost libraries
 - Suitesparse library


Here, the path to the suitesparse library is set to "/usr/include/suitesparse" in the CMakeLists. Before building makefiles, set the appropriate path according to your system!  

After compilation using cmake and make, a target called "s2d" is produced 
which can be executed by: s2d -f [path to]/square1.obj

By this the simulation is started on the reference domain Omega=[0,1]^2.

The main code of the simulation is provided in the executable 
file "pde.cpp". Here the relevant parameters are described 
in lines 18-28, where e.g. the scaling parameter sigma for the non-local repulsive term 
and the interfacial parameter epsilon can be changed 
in order to simulate certain gradient flows. 

The default configuration of the parameters is appropriate 
to the second experiment provided in the thesis Section 4.
