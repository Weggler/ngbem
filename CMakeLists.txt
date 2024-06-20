cmake_minimum_required(VERSION 3.16)

project(ngbem)
include(ngsolve_addon.cmake)

add_ngsolve_addon(_ngbem src/ngbem.cpp src/python_bem.cpp
      src/hmat.cpp src/intrules.cpp src/test_compression.cpp)

install(TARGETS _ngbem DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)
install(FILES src/__init__.py DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)

ngsolve_generate_stub_files(ngbem)
