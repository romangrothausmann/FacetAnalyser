################################################################################
# base system
################################################################################
FROM ubuntu:18.04 as system


################################################################################
# builder
################################################################################
FROM ubuntu:18.04 as builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    ca-certificates `# essential for git over https` \
    build-essential \
    cmake

### VTK
RUN git clone --depth 1 -b planeIDs4vtkHull_v8.1.2 https://gitlab.kitware.com/romangrothausmann/vtk.git

RUN mkdir -p VTK_build && \
    cd VTK_build && \
    cmake \
    	  -DCMAKE_INSTALL_PREFIX=/opt/vtk/ \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_SHARED_LIBS=ON \
	  -DBUILD_TESTING=OFF \
	  -DVTK_Group_Qt=OFF \
	  -DVTK_Group_Rendering=OFF \
	  -DVTK_RENDERING_BACKEND=None \
	  ../vtk && \
    make -j"$(nproc)" && \
    make -j"$(nproc)" install

### ITK with VTK for ITKVtkGlue
RUN git clone https://github.com/InsightSoftwareConsortium/ITK.git `# ITK-PR #765 merged into master with ITK @ dc4419da ` && \
    cd ITK && \
    git checkout dc4419daa5fa3d1e3d9ff6d8d6d76902e5d1bee9

RUN mkdir -p ITK_build && \
    cd ITK_build && \
    cmake \
    	  -DCMAKE_INSTALL_PREFIX=/opt/itk/ \
	  -DCMAKE_MODULE_PATH=/opt/vtk/lib/cmake \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_SHARED_LIBS=ON \
	  -DBUILD_TESTING=OFF \
	  -DModule_ITKVtkGlue=ON \
	  ../ITK && \
    make -j"$(nproc)" && \
    make -j"$(nproc)" install

### FacetAnalyser
COPY code/ /code/

RUN mkdir -p FacetAnalyser_build && \
    cd FacetAnalyser_build && \
    cmake \
    	  -DCMAKE_INSTALL_PREFIX=/opt/FacetAnalyser/ \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_PLUGIN=OFF \
	  -DBUILD_EXAMPLE=ON \
	  -DBUILD_TESTING=OFF \
	  ../code/ && \
    make -j"$(nproc)" && \
    mkdir -p /opt/FacetAnalyser/bin/ && cp FacetAnalyserCLI /opt/FacetAnalyser/bin/



################################################################################
# merge
################################################################################
FROM system

COPY --from=builder /opt/vtk/ /opt/vtk/
COPY --from=builder /opt/itk/ /opt/itk/
COPY --from=builder /opt/FacetAnalyser/ /opt/FacetAnalyser/

ENV LD_LIBRARY_PATH "/opt/vtk/lib/:/opt/itk/lib/:${LD_LIBRARY_PATH}"
ENV PV_PLUGIN_PATH "/opt/FacetAnalyser/"
ENV PATH "/opt/FacetAnalyser/bin/:${PATH}"
CMD ["/opt/FacetAnalyser/bin/FacetAnalyserCLI"]
