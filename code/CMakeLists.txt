
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Change PROJECT_NAME to the name of your project
PROJECT(FacetAnalyser)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(CTest)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags") # use RPATH instead of RUNPATH (currently ignored on Debian-9), -WL tells gcc to pass option to the linker

OPTION( BUILD_PLUGIN
  "Build ParaView plugin"
  ON )

OPTION( BUILD_EXAMPLE
  "Build example"
  ON )

OPTION( BUILD_TESTING
  "Build tests"
  OFF )


FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
    
IF(BUILD_PLUGIN)
FIND_PACKAGE(ParaView)
IF (${ParaView_FOUND})
  INCLUDE(${PARAVIEW_USE_FILE})
ELSE()
  SET(BUILD_PLUGIN OFF)
ENDIF()
ENDIF()

IF(BUILD_TESTING)

  ## build test programs
  ADD_EXECUTABLE(FacetAnalyserTest
    FacetAnalyser.cxx vtkGaussianSplatterExtended.cxx
    FacetAnalyserTest.cxx)
  TARGET_LINK_LIBRARIES(FacetAnalyserTest PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  
  ## run tests
  ADD_TEST(NAME basicTest01 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/ COMMAND FacetAnalyserTest 8faced-rhombic-dodecahedron_twinned_simplified.vtk)
    
ENDIF(BUILD_TESTING)

IF(BUILD_PLUGIN)

  # Use the PLUGIN macro to create a plugin.
  ADD_PARAVIEW_PLUGIN(FacetAnalyser "1.0"
    SERVER_MANAGER_SOURCES FacetAnalyser.cxx vtkGaussianSplatterExtended.cxx
    SERVER_MANAGER_XML FacetAnalyser.xml)
  
  TARGET_LINK_LIBRARIES(FacetAnalyser PRIVATE ${ITK_LIBRARIES})
  INSTALL(TARGETS FacetAnalyser LIBRARY DESTINATION . )
  
  IF(BUILD_TESTING)
    
    ADD_TEST(NAME basicPluginTest01 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/demos/ COMMAND ${ParaView_DIR}/bin/pvpython ${PROJECT_SOURCE_DIR}/scripts/pvsm2webgl.py -i 8faced-rhombic-dodecahedron_twinned_simp+def.pvsm -o ${PROJECT_SOURCE_DIR}/docs/demos/WebGL/8faced-rhombic-dodecahedron_twinned_simp+def.webgl -p ${CMAKE_CURRENT_BINARY_DIR}/libFacetAnalyser.so )
    
    ADD_TEST(NAME basicPluginTest02 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/demos/ COMMAND ${ParaView_DIR}/bin/pvpython ${PROJECT_SOURCE_DIR}/scripts/pvsm2x3d.py -i 8faced-rhombic-dodecahedron_twinned_simp+def.pvsm -o ../docs/demos/X3D/8faced-rhombic-dodecahedron_twinned_simp+def.x3d -p ${CMAKE_CURRENT_BINARY_DIR}/libFacetAnalyser.so ) # relative path as it is stored in X3D
    
    set_tests_properties(basicPluginTest01 basicPluginTest02 PROPERTIES PASS_REGULAR_EXPRESSION "Export succeded")    ## regex-check because script itself fails with: Inconsistency detected by ld.so: dl-close.c: _dl_close: Assertion `map->l_init_called' failed!
    
  ENDIF(BUILD_TESTING)
ENDIF(BUILD_PLUGIN)
  
IF(BUILD_EXAMPLE)
    
  # EXAMPLE ---------------------------------------------------------
  
  ADD_EXECUTABLE(FacetAnalyserCLI
    FacetAnalyser.cxx vtkGaussianSplatterExtended.cxx
    FacetAnalyserCLI.cxx)
  TARGET_LINK_LIBRARIES(FacetAnalyserCLI PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  
  #------------------------------------------------------------------

  IF(BUILD_TESTING)
    
    ADD_TEST( NAME basicExampleTest01 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/demos/ COMMAND FacetAnalyserCLI 8faced-rhombic-dodecahedron_twinned_010_deformed_wss.vtp 51 15 0.2 0.001 2 8faced-rhombic-dodecahedron_twinned_010_deformed_wss_FA.vtp )
    
  ENDIF(BUILD_TESTING)
ENDIF(BUILD_EXAMPLE)
