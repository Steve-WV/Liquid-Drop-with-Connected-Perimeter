cmake_minimum_required(VERSION 3.5)
project(s3d CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# My optimization flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native -DNDEBUG")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
endif(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-Ofast" COMPILER_OPT_OFAST_SUPPORTED)
if (COMPILER_OPT_OFAST_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif(COMPILER_OPT_OFAST_SUPPORTED)

find_path(CHOLMOD_INCLUDE_DIR NAMES cholmod.h PATHS "/usr/include/suitesparse" NO_DEFAULT_PATH )
if ( CHOLMOD_INCLUDE_DIR )
  list ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
endif( CHOLMOD_INCLUDE_DIR )


find_path( SUITESPARSE_LIBRARY_DIR NAMES SuiteSparse_config.h PATHS "/usr/include/suitesparse" NO_DEFAULT_PATH )
if ( SUITESPARSE_LIBRARY_DIR )
  list ( APPEND SUITESPARSE_LIBRARIES amd)
  list ( APPEND SUITESPARSE_LIBRARIES btf)
  list ( APPEND SUITESPARSE_LIBRARIES camd)
  list ( APPEND SUITESPARSE_LIBRARIES ccolamd)
  list ( APPEND SUITESPARSE_LIBRARIES cholmod)
  list ( APPEND SUITESPARSE_LIBRARIES colamd)
 # list ( APPEND SUITESPARSE_LIBRARIES csparse)
  list ( APPEND SUITESPARSE_LIBRARIES cxsparse)
  list ( APPEND SUITESPARSE_LIBRARIES klu)
  list ( APPEND SUITESPARSE_LIBRARIES spqr)
  list ( APPEND SUITESPARSE_LIBRARIES umfpack)


#find_library( SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig PATHS "/usr/include/suitesparse" NO_DEFAULT_PATH )
 find_library( SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig PATHS ${SUITESPARSE_LIBRARY_DIR})
 if ( EXISTS ${SUITESPARSE_CONFIG_LIB} )
    list ( APPEND SUITESPARSE_LIBRARIES ${SUITESPARSE_CONFIG_LIB} )
  endif ( EXISTS ${SUITESPARSE_CONFIG_LIB} )

else( SUITESPARSE_LIBRARY_DIR )  
  MESSAGE ( WARNING "SUITESPARSE_LIBRARY_DIR not found!" )
endif( SUITESPARSE_LIBRARY_DIR )  

if (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set(SUITESPARSE_FOUND TRUE)
else (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set( SUITESPARSE_FOUND FALSE )
endif (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

include_directories(${SUITESPARSE_INCLUDE_DIRS})
link_directories(${SUITESPARSE_LIBRARY_DIR})
link_libraries(${SUITESPARSE_LIBRARIES})



# Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
set(EIGEN_USE_BLAS=1)
set(EIGEN_USE_LAPACKE=1)
set(EIGEN_USE_UMFPACK=1)

# VTK
find_package(VTK)
if (VTK_FOUND)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_VTK")
   include(${VTK_USE_FILE})
   link_libraries(${VTK_LIBRARIES})
endif(VTK_FOUND)

# BLAS
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})
link_libraries(${BLAS_LIBRARIES})

# LAPACKE
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})
link_libraries(${LAPACK_LIBRARIES})


# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# OpenMP
find_package(OpenMP)
if (OpenMP_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)


add_executable(s2d mesh.cpp connected.cpp pde.cpp export.cpp main.cpp)
