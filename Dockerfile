FROM registry.gitlab.com/romangrothausmann/dockerfiles/vtk-rhg/itk:vtk-7.1.1_planeIDs4vtkHull_itk-4.12.2 as builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    ca-certificates `# essential for git over https` \
    cmake \
    build-essential

RUN git clone http://github.com/romangrothausmann/FacetAnalyser.git

RUN mkdir -p FacetAnalyser_build && \
    cd FacetAnalyser_build && \
    cmake \
    	  -DCMAKE_INSTALL_PREFIX=/opt/FacetAnalyser/ \
	  -DCMAKE_MODULE_PATH=/opt/vtk/lib/cmake \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_PLUGIN=OFF \
	  -DBUILD_EXAMPLE=ON \
	  -DBUILD_TESTING=OFF \
	  ../FacetAnalyser/code/ && \
    make -j"$(nproc)" && \
    make -j"$(nproc)" install


FROM ubuntu:16.04

COPY --from=builder /opt/vtk/ /opt/vtk/
COPY --from=builder /opt/itk/ /opt/itk/
COPY --from=builder /opt/FacetAnalyser/ /opt/FacetAnalyser/
