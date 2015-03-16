# FacetAnalyser
ParaView plugin for automated facet detection and measurement of interplanar angles of tomographic objects

See http://www.midasjournal.org/browse/publication/951 or http://hdl.handle.net/10380/3510 for more details.



## Installation


Configure ParaView with cmake as follows:

`BUILD_SHARED_LIBS  ON`  

Compile ParaView with make and optionally install it. The following ITK compilation does not need ParaView to be installed. After compilation change into the directory /paraview-build-dir/lib/ and create symbolic links without the ParaView suffix. In a BASH fore example with:

`for i in *-pv*.so; do ln -s $i ${i%-pv*}.so; done`


Configure ITK with cmake as follows, set /paraview-build-dir/ to the build directory used for building ParaView:


`Module_ITKVtkGlue  ON`  
`Module_ITKReview   ON`   
`VTK_DIR            /paraview-build-dir/VTK/`  
`CMAKE_CXX_FLAGS    -L/paraview-build-dir/lib/`  

The two additionally enabled ITK modules are needed for the connection of VTK with ITK and for the watershed filters. It is essential, that VTK_DIR is set to the build directory containing VTK shipped with ParaView.

