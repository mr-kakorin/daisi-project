#ifndef BOUNDARYCONTAINER2D_H
#define BOUNDARYCONTAINER2D_H

#include <common_tools/constants.h>

#include "Dmath.h"
#include "Geom.h"
#include "VTKIncludeSolver.h"
#include <cmath>

template <class PointType>
class BoundaryContainer2dpolar
{
    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& xmin;
        ar& xmax;
        ar& ymin;
        ar& ymax;
        ar& EdgesData;
        ar& VertexData;
        ar& ContainerSize;
        ar& epsilon;
        ar& flagInit;
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& xmin;
        ar& xmax;
        ar& ymin;
        ar& ymax;
        ar& EdgesData;
        ar& VertexData;
        ar& ContainerSize;
        ar& epsilon;
        ar& flagInit;
        epsilon = 1e-12 * std::max((xmax - xmin), (ymax - ymin));
        if (flagInit != 0)
            ConvertBoundary2VTKUnstructuredGrid();
    };
    BoundaryContainer2d<PointType>& operator=(const BoundaryContainer2d<PointType>& right);

    PointType epsilon;

    PointType xmin;
    PointType xmax;
    PointType ymin;
    PointType ymax;

    std::vector<DGeo::Edge<PointType>>  EdgesData;
    std::vector<DGeo::Point<PointType>> VertexData;

    vtkSmartPointer<vtkUnstructuredGrid> VTKgrid;

    int ContainerSize;
    int flagInit;

    void RemoveFineEdges();
    void Merge(BoundaryContainer2d<PointType> obj);
    BoundaryContainer2d(std::string InputFile);
    DGeo::Point<PointType> findPoint(int n);
    BoundaryContainer2d();
    void GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr);
    void GetPlotRotate(std::vector<std::vector<float>>& Xarr,
                       std::vector<std::vector<float>>& Yarr);
    vtkSmartPointer<vtkUnstructuredGrid> GetBoundaryVTKUnstructuredGrid();
    void                                 WriteBoundary2VTK(std::string InputFileName);
    int                                  NumberIntersections(DGeo::Edge<PointType>);
    bool IsIntersection(DGeo::Edge<PointType>, DGeo::Point<PointType>& intersectionPoint,
                        int&) const;
    void FindAllIntersections(DGeo::Edge<PointType>                edge,
                              std::vector<DGeo::Point<PointType>>& intersectionPoint,
                              std::vector<DGeo::Edge<PointType>>&  intersectionEdgesTmp) const;
    int  InBoundary(DGeo::Point<PointType>, PointType h);
    void MaxVertex();
    int  OnBoundary(DGeo::Point<PointType>);
    void ConvertBoundary2VTKUnstructuredGrid();
    bool IsIntersectionLight(DGeo::Edge<PointType> edge, DGeo::Point<PointType>& intersectionPoint,
                             int& intersectionEdgeNumber) const;
};

#endif