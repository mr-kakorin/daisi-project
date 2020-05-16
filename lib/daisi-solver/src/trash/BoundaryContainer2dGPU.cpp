#include "BoundaryContainer2dGPU.h"
template BoundaryContainer2dGPU<float>;
template BoundaryContainer2dGPU<double>;
template <class PointType>
BoundaryContainer2dGPU<PointType>::BoundaryContainer2dGPU(BoundaryContainer2d<PointType>& obj)
{
    xmin = obj.xmin;
    xmax = obj.xmax;
    ymin = obj.ymin;
    ymax = obj.ymax;

    EdgesData  = obj.EdgesData;
    VertexData = obj.VertexData;

    ContainerSize = obj.ContainerSize;
    flagInit      = obj.flagInit;
};