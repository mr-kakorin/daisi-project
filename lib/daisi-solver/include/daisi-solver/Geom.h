#ifndef GEOM_H
#define GEOM_H
#include "VTKIncludeSolver.h"
#include "boost/serialization/vector.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <vector>

namespace DGeo
{

template <class PointType> class Point
{
  friend class boost::serialization::access;

  public:
    PointType x;
    PointType y;
    PointType z;
    int       Number;
    bool      isOut;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& x;
        ar& y;
        ar& z;
        ar& Number;
        ar& isOut;
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& x;
        ar& y;
        ar& z;
        ar& Number;
        ar& isOut;
    }
    bool isEqual(Point<PointType> point1, PointType epsilon) const;
    PointType Radius();
    PointType Dist2Point(const Point<PointType>& point1) const;
    Point<PointType> rotate(PointType alpha, Point<PointType> point);
    Point<PointType>& operator=(const Point<PointType>& right);
    Point<PointType> operator+(const Point<PointType>& right) const;
    Point<PointType> operator-(const Point<PointType>& right) const;
    Point();
    Point(Point<PointType>, PointType, PointType, PointType);
};

template <typename PointType>
struct CurvePoints_t
{
    std::vector<PointType> x;
    std::vector<PointType> y;
    std::vector<PointType> z;
};

template <typename PointType>
CurvePoints_t<PointType> calc_grad2d(std::vector<DGeo::Point<PointType>> const& in_vec) {
    CurvePoints_t<PointType> result;
    result.x.resize( in_vec.size() );
    result.y.resize( in_vec.size() );
    double h = 1e-8;
    result.x[0] = (in_vec[1].x - in_vec[0].x) / h;
    result.x[result.x.size() - 1] = (in_vec[result.x.size() - 1].x - in_vec[result.x.size() - 2].x) / h;
    result.y[0] = (in_vec[1].y - in_vec[0].y) / h;
    result.y[result.y.size() - 1] = (in_vec[result.y.size() - 1].y - in_vec[result.y.size() - 2].y) / h;
    std::size_t end = in_vec.size();
    for (int i = 1; i < end; ++i) {
        result.x[i] = (in_vec[i + 1].x - in_vec[i - 1].x) / (2 * h);
        result.y[i] = (in_vec[i + 1].y - in_vec[i - 1].y) / (2 * h);
    }
    return result;
}

template <class PointType>
inline PointType Pfabs(PointType val)
{
    if (val.x < 0)
        val.x = -val.x;
    if (val.y < 0)
        val.y = -val.y;
    if (val.z < 0)
        val.z = -val.z;
    return val;
}

template <class PointType>
inline int PointCmp(Point<PointType> p1, Point<PointType> p2, int key)
{
    if (key == 0)
    {
        if ((p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z))
            return 0;
        if (p1.y < p2.y)
            return -1;
        if (p1.y == p2.y)
        {
            if (p1.x < p2.x)
                return -1;
        }
    }
    return 1;
}

template <class PointType>
class Edge
{
  public:
    Point<PointType> point1;
    Point<PointType> point2;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& point2;
        ar& point1;
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& point2;
        ar& point1;
    }
    Point<PointType> GetNormalPoint1(PointType h) const;
    Point<PointType> GetNormalPoint2(PointType h) const;

    Edge<PointType> TranslateByNormal1(PointType h) const;
    Edge<PointType> TranslateByNormal2(PointType h) const;

    vtkSmartPointer<vtkUnstructuredGrid> GetVTKGrid();
    PointType                            alpha() const;
    Point<PointType>                     Middle() const;
    std::vector<Edge<PointType>> resize(int size);
    Edge<PointType>& operator=(const Edge<PointType>& right);
    PointType length() const;
    bool IsParallel(const DGeo::Edge<PointType>& edge2, double epsilon) const;
    bool IsEqual(const DGeo::Edge<PointType>& edge2, double epsilon) const;
    bool IsParallelFixEps(const DGeo::Edge<PointType>& edge2, double epsilon) const;
    bool IsIntersectionEdge(const DGeo::Edge<PointType>& edge2, double epsilon) const;
    bool IsIntersectionEdgeMesh(const DGeo::Edge<PointType>& edge2, double epsilon) const;

    bool IsIntersectionEdge(const DGeo::Edge<PointType>& edge2, double epsilon,
                            DGeo::Point<PointType>* out_intersection) const;
    bool IsIntersectionEdgeLight(const DGeo::Edge<PointType>& edge2, double epsilon,
                                 DGeo::Point<PointType>* out_intersection) const;
    bool IsIntersectionEdgePoints(const DGeo::Edge<PointType>& edge2, double epsilon,
                                  DGeo::Point<PointType>* out_intersection) const;
    bool IsOnOneSide(const DGeo::Point<PointType>& pp1, const DGeo::Point<PointType>& pp2, double epsilonLoc) const;
    int IsOnEdge(const DGeo::Point<PointType>& point, double epsilon) const;
    bool IsOnLine(const DGeo::Point<PointType>& point, double epsilon) const;
    void SwapPoints();
};
}
#endif