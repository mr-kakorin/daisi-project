#ifndef BOUNDARYCONTAINER3D_H
#define BOUNDARYCONTAINER3D_H
#include "VTKIncludeSolver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
}
template <class PointType>
class BoundaryContainer3d
{
    friend class boost::serialization::access;

  public:
    std::vector<DGeo::Edge<PointType>> EdgesData;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    void RemoveFineEdges(float k);
    BoundaryContainer3d();
    BoundaryContainer3d(std::string InputFile, std::string& errorMsg);
    BoundaryContainer3d(std::vector<std::vector<double>> vertex, bool IsTube,
                        double roundingRadius);
    BoundaryContainer3d(std::vector<std::vector<double>> vertex, bool IsFirst, bool IsLast,
                        double roundingRadius);
    void GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr);
    void GetPlotRotate(std::vector<std::vector<float>>& Xarr, std::vector<std::vector<float>>& Yarr,
                       int angles);
    void GetPlotRotate1(std::vector<std::vector<float>>& Xarr,
                        std::vector<std::vector<float>>& Yarr, int angles);
    void ConvertBoundary2VTKUnstructuredGrid();
    void Merge(const std::shared_ptr<BoundaryContainer3d<PointType>> obj);
    vtkSmartPointer<vtkUnstructuredGrid> GetBoundaryVTKUnstructuredGrid();
};

#endif