#ifndef BOUNDARYCONTAINER2DGPU_H
#define BOUNDARYCONTAINER2DGPU_H
#include "BoundaryContainer2d.h"
#include "Dmath.h"
#include "GPUTypes.h"
#include "Geom.h"
#include "GeomGPU.h"
#include "VTKIncludeSolver.h"
#include <common_tools/constants.h>
#include <cmath>
#include <cuda_runtime.h>
template <class PointType> class BoundaryContainer2dGPU
{
  public:
    PointType epsilon;

    PointType xmin;
    PointType xmax;
    PointType ymin;
    PointType ymax;

    GPUTypes::vector<DGeoGPU::Edge<PointType>>  EdgesData;
    GPUTypes::vector<DGeoGPU::Point<PointType>> VertexData;

    int ContainerSize;
    int flagInit;

    BoundaryContainer2dGPU(){

    };
    BoundaryContainer2dGPU(const std::shared_ptr<BoundaryContainer2d<PointType>> obj);
    /*BoundaryContainer2dGPU(BoundaryContainer2d<PointType> boundaryCPU);
    bool IsIntersection(DGeoGPU::Edge<PointType>);
    bool IsIntersection(DGeoGPU::Edge<PointType>, DGeoGPU::Point<PointType>& intersectionPoint, int&) const;
    int InBoundary(DGeoGPU::Point<PointType>, PointType h);
    void MaxVertex();
    int OnBoundary(DGeoGPU::Point<PointType>);*/
};

#endif