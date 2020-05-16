#ifndef BOUNDARYCONTAINER2D_H
#define BOUNDARYCONTAINER2D_H

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
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;

template <class PointType>
class BoundaryContainer2d
{
    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    BoundaryContainer2d<PointType>& operator=(const BoundaryContainer2d<PointType>& right);

    PointType              xmin;
    PointType              xmax;
    PointType              ymin;
    PointType              ymax;
    DGeo::Point<PointType> PointTmp;
    DGeo::Edge<PointType>  EdgeTmp;

    std::vector<DGeo::Edge<PointType>>  EdgesData;
    std::vector<DGeo::Point<PointType>> VertexData;

    vtkSmartPointer<vtkUnstructuredGrid> VTKgrid;

    int  ContainerSize;
    int  flagInit;
    void ttf(std::string& InputFile, std::string& errorMsg);

    void sort();
    void Reverse();
    void RemoveFineEdges(float k);
    void Merge(const std::shared_ptr<BoundaryContainer2d<PointType>> obj);
    // int MergeWithSort(std::shared_ptr <BoundaryContainer2d<PointType>> obj);
    BoundaryContainer2d(std::string& InputFile, std::string& errorMsg);
    BoundaryContainer2d(std::vector<std::vector<double>> vertex, bool IsTube,
                        double roundingRadius); // one period
    BoundaryContainer2d(std::vector<std::vector<double>> vertex, bool IsFirst, bool IsLast,
                        double roundingRadius); // many periods
    void AngleRounding(std::vector<DGeo::Point<PointType>>& vertex, double roundingRadius,
                       int& num); //����������� �������� �� 3� �������� ����
    DGeo::Point<PointType> findPoint(int n);
    BoundaryContainer2d();
    void GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr);
    void GetPlotRotate(std::vector<std::vector<float>>& Xarr, std::vector<std::vector<float>>& Yarr,
                       int angles);
    void GetPlotRotate1(std::vector<std::vector<float>>& Xarr,
                        std::vector<std::vector<float>>& Yarr, int angles);
    void AddEdge(DGeo::Edge<PointType> edg);
    void AddVertex(DGeo::Point<PointType> ver);
    vtkSmartPointer<vtkUnstructuredGrid> GetBoundaryVTKUnstructuredGrid();
    void WriteBoundary2VTK(std::string InputFileName);
    int NumberIntersections(double epsilon, const DGeo::Edge<PointType>&);
    bool IsIntersection(double                  epsilon, DGeo::Edge<PointType>,
                        DGeo::Point<PointType>& intersectionPoint, int&) const;
    void FindAllIntersections(double epsilon, DGeo::Edge<PointType> edge,
                              std::vector<DGeo::Point<PointType>>& intersectionPoint,
                              std::vector<DGeo::Edge<PointType>>&  intersectionEdgesTmp) const;
    int InBoundary(double epsilon, DGeo::Point<PointType>, PointType h);
    void MaxVertex();
    int OnBoundary(double epsilon, DGeo::Point<PointType>);
    void ConvertBoundary2VTKUnstructuredGrid();
    bool IsIntersectionLight(double epsilon, DGeo::Edge<PointType> edge,
                             DGeo::Point<PointType>& intersectionPoint,
                             int&                    intersectionEdgeNumber) const;
};
#endif