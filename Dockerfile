FROM registry.gitlab.com/romangrothausmann/dockerfiles/vtk-rhg/itk:vtk-7.1.1_planeIDs4vtkHull_itk-4.12.2 as builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    build-essential \
    libsm-dev libx11-dev libxt-dev libxext-dev `# needed for ITKVtkGlue`

COPY code/ /code/

RUN mkdir -p FacetAnalyser_build && \
    cd FacetAnalyser_build && \
    cmake \
    	  -DCMAKE_INSTALL_PREFIX=/opt/FacetAnalyser/ \
	  -DCMAKE_PREFIX_PATH='/opt/vtk/lib/cmake/;/opt/itk/lib/cmake/' \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_PLUGIN=OFF \
	  -DBUILD_EXAMPLE=ON \
	  -DBUILD_TESTING=OFF \
	  ../code/ && \
    make -j"$(nproc)" && \
    mkdir -p /opt/FacetAnalyser/bin/ && cp FacetAnalyserCLI /opt/FacetAnalyser/bin/


FROM ubuntu:16.04

COPY --from=builder /opt/vtk/ /opt/vtk/
COPY --from=builder /opt/itk/ /opt/itk/
COPY --from=builder /opt/FacetAnalyser/ /opt/FacetAnalyser/

ENV PATH="/opt/FacetAnalyser/bin/:${PATH}"
CMD ["/opt/FacetAnalyser/bin/FacetAnalyserCLI"]
