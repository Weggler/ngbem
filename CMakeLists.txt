cmake_minimum_required(VERSION 3.28)

include(FetchContent)

project(ngbem)

option( USE_KiFMM "Use KiFMM" OFF )
option( USE_FMM3D "Use FMM3D" OFF )

include(ngsolve_addon.cmake)

add_ngsolve_addon(_ngbem src/ngbem.cpp src/python_bem.cpp
      src/hmat.cpp src/intrules.cpp src/test_compression.cpp)

if(USE_KiFMM)
  target_compile_definitions(_ngbem PRIVATE USE_KiFMM)
  FetchContent_Declare(
    kifmm
    GIT_REPOSITORY https://github.com/bempp/kifmm.git
    GIT_TAG "enh/c-abi"
    EXCLUDE_FROM_ALL
  )
  FetchContent_Populate(kifmm)
  add_custom_target(build_kifmm ALL
    COMMAND cargo build -r WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/kifmm-src
    BYPRODUCTS ${FETCHCONTENT_BASE_DIR}/kifmm-src/targets/kifmm_rs.h
    USES_TERMINAL
  )

  target_include_directories(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/kifmm-src/kifmm/include)
  target_link_libraries(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/kifmm-src/target/release/libkifmm.so)
  install(FILES ${FETCHCONTENT_BASE_DIR}/kifmm-src/target/release/libkifmm.so DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)
  add_dependencies(_ngbem build_kifmm)

  set_target_properties(_ngbem PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${ADDON_INSTALL_DIR_PYTHON}/ngbem/"
  )
  
else(USE_KiFMM)
  if(USE_FMM3D)
    target_compile_definitions(_ngbem PRIVATE USE_FMM3D)
    FetchContent_Declare(
      fmm3d
      GIT_REPOSITORY https://github.com/flatironinstitute/FMM3D.git
      GIT_TAG master
      EXCLUDE_FROM_ALL
    )
    FetchContent_Populate(fmm3d)

    set(fmm_lib_path ${fmm3d_SOURCE_DIR}/lib-static/libfmm3d.a)

    add_custom_command(OUTPUT ${fmm_lib_path}
      COMMAND make lib -j 9 WORKING_DIRECTORY ${fmm3d_SOURCE_DIR}
      USES_TERMINAL
    )
    add_custom_target(build_fmm DEPENDS ${fmm_lib_path})

    add_library(fmm STATIC IMPORTED GLOBAL)
    set_target_properties(fmm PROPERTIES IMPORTED_LOCATION ${fmm_lib_path})
    add_dependencies(fmm build_fmm)

    find_library(fortran_lib NAMES gfortran)
    find_package(OpenMP)

    target_include_directories(_ngbem PRIVATE ${fmm3d_SOURCE_DIR}/fmm3d-src/c)
    target_link_libraries(_ngbem PRIVATE fmm ${fortran_lib} ${OpenMP_omp_LIBRARY})
  endif(USE_FMM3D)
endif(USE_KiFMM)

install(TARGETS _ngbem DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)
install(FILES src/__init__.py DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)

ngsolve_generate_stub_files(ngbem)
