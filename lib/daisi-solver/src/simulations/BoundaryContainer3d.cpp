#include <common_tools/constants.h>

#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "Geom.h"
#include <cmath>

template class BoundaryContainer3d<float>;
template class BoundaryContainer3d<double>;

template <class PointType>
void BoundaryContainer3d<PointType>::RemoveFineEdges(float k){

};

template <class PointType>
BoundaryContainer3d<PointType>::BoundaryContainer3d(){

};
template <class PointType>
BoundaryContainer3d<PointType>::BoundaryContainer3d(std::string InputFile, std::string& errorMsg){

};
template <class PointType>
BoundaryContainer3d<PointType>::BoundaryContainer3d(std::vector<std::vector<double>> vertex,
                                                    bool IsTube, double roundingRadius){

};
template <class PointType>
BoundaryContainer3d<PointType>::BoundaryContainer3d(std::vector<std::vector<double>> vertex,
                                                    bool IsFirst, bool IsLast,
                                                    double roundingRadius){

};

template <class PointType>
void BoundaryContainer3d<PointType>::GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr){

};
template <class PointType>
void BoundaryContainer3d<PointType>::GetPlotRotate(std::vector<std::vector<float>>& Xarr,
                                                   std::vector<std::vector<float>>& Yarr,
                                                   int                              angles){};
template <class PointType>
void BoundaryContainer3d<PointType>::GetPlotRotate1(std::vector<std::vector<float>>& Xarr,
                                                    std::vector<std::vector<float>>& Yarr,
                                                    int                              angles){};
template <class PointType>
void BoundaryContainer3d<PointType>::ConvertBoundary2VTKUnstructuredGrid(){

};
template <class PointType>
void BoundaryContainer3d<PointType>::Merge(
    const std::shared_ptr<BoundaryContainer3d<PointType>> obj){

};

template <class PointType>
vtkSmartPointer<vtkUnstructuredGrid>
BoundaryContainer3d<PointType>::GetBoundaryVTKUnstructuredGrid()
{
    return NULL;
}

template void BoundaryContainer3d<float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void BoundaryContainer3d<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void BoundaryContainer3d<double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void BoundaryContainer3d<float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void BoundaryContainer3d<PointType>::save(Archive& ar, const unsigned int) const {};
template <class PointType>
template <class Archive>
void BoundaryContainer3d<PointType>::load(Archive& ar, const unsigned int){};