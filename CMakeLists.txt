cmake_minimum_required(VERSION 3.16)

###############################################################################
# Edit here
###############################################################################

# make sure the module_name matches the first argument of PYBIND11_MODULE
set(module_name libbem)
set(source_files src/ngbem.cpp src/python_bem.cpp
      src/hmat.cpp src/intrules.cpp src/test_compression.cpp)

###############################################################################
# STOP editing here
###############################################################################

# Don't change this
project(${module_name})
include(ngsolve_addon.cmake)
