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

  target_compile_definitions(_ngbem PRIVATE USE_FMM3D)
  target_include_directories(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/kifmm-src/target)
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

    # TODO: Add BYPRODUCTS to track if is up to date
    add_custom_target(build_fmm3d ALL
      COMMAND make lib -j 9 WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/fmm3d-src
      USES_TERMINAL
    )

    add_dependencies(_ngbem build_fmm3d)
    target_include_directories(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/fmm3d-src/c)
    target_link_libraries(_ngbem PRIVATE ${FETCHCONTENT_BASE_DIR}/fmm3d-src/lib-static/libfmm3d.a)
  endif(USE_FMM3D)
endif(USE_KiFMM)

install(TARGETS _ngbem DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)
install(FILES src/__init__.py DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/ngbem/)

ngsolve_generate_stub_files(ngbem)
