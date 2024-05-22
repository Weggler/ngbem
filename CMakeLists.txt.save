cmake_minimum_required(VERSION 3.5)
project(ngbem)
 
find_package(NGSolve CONFIG REQUIRED
  HINTS /usr $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

add_ngsolve_python_module(libbem src/ngbem.cpp src/python_bem.cpp
      src/hmat.cpp src/intrules.cpp src/test_compression.cpp)



# check if CMAKE_INSTALL_PREFIX is set by user, otherwise install to netgen directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NETGEN_DIR} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)

message("The python package will be installed to: ${CMAKE_INSTALL_PREFIX}/${install_dir}")
install(TARGETS libbem DESTINATION ${install_dir})
