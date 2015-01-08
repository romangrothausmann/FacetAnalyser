
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Change PROJECT_NAME to the name of your project
PROJECT(FacetAnalyser)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(CTest)


FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${USE_ITK_FILE})

INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})


SET (PARAVIEW_ADDITIONAL_LIBRARIES
   ${PARAVIEW_ADDITIONAL_LIBRARIES}
   ${ITK_LIBRARIES}
)


IF(BUILD_TESTING)

## build test programs
SET(CurrentExe "TestFacetAnalyser")
ADD_EXECUTABLE(${CurrentExe}
  FacetAnalyser.cxx 
  TestFacetAnalyser.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${ITK_LIBRARIES})

## run tests
ADD_TEST(NAME basicTest01 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/ COMMAND TestFacetAnalyser particle_0195.vtp 101 10 0.2 0.002 particle_0195_FA0.vtp particle_0195_FA1.vtp particle_0195_FA2.vtp)

ENDIF(BUILD_TESTING)


OPTION( BUILD_PLUGIN
        "Build ParaView plugin"
        ON )

IF(BUILD_PLUGIN)

# Use the PLUGIN macro to create a plugin.
ADD_PARAVIEW_PLUGIN(FacetAnalyser "1.0"
   SERVER_MANAGER_SOURCES FacetAnalyser.cxx vtkGaussianSplatterExtended.cxx
   SERVER_MANAGER_XML FacetAnalyser.xml)

TARGET_LINK_LIBRARIES(FacetAnalyser ${ITK_LIBRARIES})

ENDIF(BUILD_PLUGIN)
