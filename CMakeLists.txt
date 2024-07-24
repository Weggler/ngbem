cmake_minimum_required(VERSION 3.28)

include(FetchContent)

project(ngbem)

option( USE_FMM3D "Use FMM3D" OFF )

include(ngsolve_addon.cmake)

add_ngsolve_addon(_ngbem src/ngbem.cpp src/python_bem.cpp
      src/hmat.cpp src/intrules.cpp src/test_compression.cpp)

if(USE_FMM3D)
  target_compile_definitions(_ngbem PRIVATE USE_FMM3D)
  FetchContent_Declare(
    fmm3d
    GIT_REPOSITORY https://github.com/flatironinstitute/FMM3D.git
    GIT_TAG master
    EXCLUDE_FROM_ALL
  )
  FetchContent_Populate(fmm3d)

  # TODO: Add BYPRODUCTS to track if is up to date
  add_custom_target(build_fmm3d ALL
    COMMAND make lib -j 9 WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/fmm3d-src
    USES_TERMINAL
  )

  add_dependencies(_ngbem build_fmm3d)
  target_include_directories(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/fmm3d-src/c)
  target_link_libraries(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/fmm3d-src/lib-static/libfmm3d.a)
endif(USE_FMM3D)

install(TARGETS _ngbem DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)
install(FILES src/__init__.py DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)

ngsolve_generate_stub_files(ngbem)
