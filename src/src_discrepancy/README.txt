This parallel OpenMP-based code has been developed by Quientin FERRE, based on Adrien Pilleboue's non-parallel implementation.

Presently (April 2016), it run under Linuxusing OpenMP
On Mac, OpenMP is disabled because of non-compatibility of OpenMP and boost

To enabla/disablee OpenMP, make 2 things:

1. in CMakeLists.txt: comment/uncomment find_package(OpenMP REQUIRED) and ${OpenMP_CXX_FLAGS}
2. it the very beginning)n src/parallelism.cpp: comment/uncomment #define  OMP_ENABLED 

Copy CMakeLists_LINUX.txt or CMakeLists_MAC.txt -> CMakeLists.txt
Copy src/parallelism_LINUX.cpp or src/parallelism_MAC.cpp -> src/parallelism.cpp