#ifndef GEOMGPU_H
#define GEOMGPU_H
#include "Geom.h"
#include "VTKIncludeSolver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>
namespace DGeoGPU
{

template <class PointType> class Point
{
  public:
    PointType x;
    PointType y;
    PointType z;
    int       Number;
    Point(DGeo::Point<PointType> obj)
    {
        x      = obj.x;
        y      = obj.y;
        z      = obj.z;
        Number = obj.Number;
    };
    PointType Radius()
    {
        return sqrt(x * x + y * y + z * z);
    };
    Point<PointType> rotate(PointType alpha, Point<PointType> point)
    {
        Point<PointType> result = *this;
        result                  = result - point;

        PointType newx = std::cos(alpha) * result.x - std::sin(alpha) * result.y;
        PointType newy = std::sin(alpha) * result.x + std::cos(alpha) * result.y;

        result.x = newx;
        result.y = newy;

        result = result + point;
        return result;
    };
    Point<PointType>& operator=(const Point<PointType>& right)
    {
        if (this == &right)
        {
            return *this;
        }
        x      = right.x;
        y      = right.y;
        z      = right.z;
        Number = right.Number;
        return *this;
    };
    Point<PointType> operator+(const Point<PointType>& right) const
    {
        Point<PointType> result = *this;
        result.x += right.x;
        result.y += right.y;
        result.z += right.z;
        return result;
    };
    Point<PointType> operator+(const Point<PointType>& right)
    {
        Point<PointType> result = *this;
        result.x += right.x;
        result.y += right.y;
        result.z += right.z;
        return result;
    };
    Point<PointType> operator-(const Point<PointType>& right) const
    {
        Point<PointType> result = *this;
        result.x -= right.x;
        result.y -= right.y;
        result.z -= right.z;
        return result;
    };
    Point<PointType> operator-(const Point<PointType>& right)
    {
        Point<PointType> result = *this;
        result.x -= right.x;
        result.y -= right.y;
        result.z -= right.z;
        return result;
    };
    Point();
    Point<PointType>* AllocateOnGPU()
    {
        Point<PointType>* result;
        cudaMalloc((void**)&result, sizeof(Point<PointType>));
        cudaMemcpy(result, this, sizeof(Point<PointType>), cudaMemcpyHostToDevice);
        return result;
    };
    Point(Point<PointType>, PointType, PointType, PointType);
};
template <class PointType> inline PointType Pfabs(PointType val)
{
    if (val.x < 0)
        val.x = -val.x;
    if (val.y < 0)
        val.y = -val.y;
    if (val.z < 0)
        val.z = -val.z;
    return val;
};
template <class PointType>
Point<PointType>::Point(){

};
template <class PointType> Point<PointType>::Point(Point<PointType> Old, PointType hx, PointType hy, PointType hz)
{
    *this = Old;
    x += hx;
    y += hy;
    z += hz;
};
template <class PointType> int PointCmp(Point<PointType> p1, Point<PointType> p2, int key)
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
        };
    };
    return 1;
};

template <class PointType> class Edge
{
  public:
    Point<PointType> point1;
    Point<PointType> point2;
    Edge(){

    };
    Edge(DGeo::Edge<PointType> obj)
    {
        point1 = obj.point1;
        point2 = obj.point2;
    };

    PointType alpha() const
    {
        if (point1.x == point2.x)
            return commtools::PI() / 2;
        return atan((point1.y - point2.y) / ((point1.x - point2.x)));
    }
    Point<PointType> Middle()
    {
        Point<PointType> p;
        p.x = (point1.x + point2.x) / 2;
        p.y = (point1.y + point2.y) / 2;
        p.z = (point1.z + point2.z) / 2;
        return p;
    };
    std::vector<Edge<PointType>> resize(int size)
    {
        std::vector<Edge<PointType>> result;
        Point<PointType>             pointTmp = point1;
        PointType                    dx       = (point2.x - point1.x) / size;
        PointType                    dy       = (point2.y - point1.y) / size;
        PointType                    dz       = (point2.z - point1.z) / size;
        for (int i = 0; i < size; i++)
        {
            Edge<PointType>* newEdge = new Edge<PointType>();
            newEdge->point1          = pointTmp;

            newEdge->point2.x = newEdge->point1.x + dx;
            newEdge->point2.y = newEdge->point1.y + dy;
            newEdge->point2.z = newEdge->point1.z + dz;

            pointTmp = newEdge->point2;
            result.push_back(*newEdge);
        };
        result.back().point2.x = point2.x;
        result.back().point2.y = point2.y;
        result.back().point2.z = point2.z;
        return result;
    };
    Edge<PointType>& operator=(const Edge<PointType>& right)
    {
        if (this == &right)
        {
            return *this;
        }
        point1 = right.point1;
        point2 = right.point2;
        return *this;
    };
    PointType length()
    {
        return sqrt((point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y) +
                    (point1.z - point2.z) * (point1.z - point2.z));
    };
    template <class PointType> bool IsParallel(DGeo::Edge<PointType> edge2, PointType epsilon)
    {
        if (std::abs((edge2.point1.x - point1.x) * (point2.y - point1.y) -
                (edge2.point1.y - point1.y) * (point2.x - point1.x)) < epsilon &&
            std::abs((edge2.point2.x - point1.x) * (point2.y - point1.y) -
                (edge2.point2.y - point1.y) * (point2.x - point1.x)) < epsilon)
            return true;
        return false;
    };
    template <class PointType> bool IsIntersectionEdge(DGeo::Edge<PointType> edge2, PointType epsilon)
    {
        DGeo::Edge<PointType>  edge1  = *this;
        DGeo::Point<PointType> start1 = edge1.point1;
        DGeo::Point<PointType> start2 = edge2.point1;
        DGeo::Point<PointType> end1   = edge1.point2;
        DGeo::Point<PointType> end2   = edge2.point2;

        if (edge1.IsParallel(edge2, epsilon))
            return false;

        if (edge1.IsOnEdge(edge2.point2, epsilon) < 2 || edge1.IsOnEdge(edge2.point1, epsilon) < 2 ||
            edge2.IsOnEdge(edge1.point1, epsilon) < 2 || edge2.IsOnEdge(edge1.point2, epsilon) < 2)
            return true;

        DGeo::Point<PointType> dir1 = Pfabs(edge1.point2 - edge1.point1);
        DGeo::Point<PointType> dir2 = Pfabs(edge2.point2 - edge2.point1);

        //������� ��������� ������ ���������� ����� �������
        PointType a1 = edge1.point1.y - edge1.point2.y;
        PointType b1 = edge1.point2.x - edge1.point1.x;
        PointType d1 = edge1.point1.x * edge1.point2.y - edge1.point2.x * edge1.point1.y;

        PointType a2 = edge2.point1.y - edge2.point2.y;
        PointType b2 = edge2.point2.x - edge2.point1.x;
        PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

        //����������� ����� ��������, ��� ��������� � ����� ������������� ���
        PointType line1_start = a2 * edge1.point1.x + b2 * edge1.point1.y + d2;
        PointType line1_end   = a2 * edge1.point2.x + b2 * edge1.point2.y + d2;

        PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
        PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

        //���� ����� ������ ������� ����� ���� ����, ������ �� � ����� ������������� � ����������� ���.
        if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
            return true;

        //	PointType u = seg1_line2_start / (seg1_line2_start - seg1_line2_end);
        //	*out_intersection = start1 + u*dir1;

        return false;
    }
    bool IsIntersectionEdge(DGeo::Edge<PointType> edge2, PointType epsilon, DGeo::Point<PointType>* out_intersection)
    {
        DGeo::Edge<PointType>  edge1  = *this;
        DGeo::Point<PointType> start1 = edge1.point1;
        DGeo::Point<PointType> start2 = edge2.point1;
        DGeo::Point<PointType> end1   = edge1.point2;
        DGeo::Point<PointType> end2   = edge2.point2;

        if (edge1.IsParallel(edge2, epsilon))
            return false;

        if (edge1.IsOnEdge(edge2.point2, epsilon) < 2)
            return true;

        DGeo::Point<PointType> dir1  = Pfabs(edge1.point2 - edge1.point1);
        DGeo::Point<PointType> dir2  = Pfabs(edge2.point2 - edge2.point1);
        DGeo::Point<PointType> dir11 = edge1.point2 - edge1.point1;

        PointType a1 = edge1.point1.y - edge1.point2.y;
        PointType b1 = edge1.point2.x - edge1.point1.x;
        PointType d1 = edge1.point1.x * edge1.point2.y - edge1.point2.x * edge1.point1.y;

        PointType a2 = edge2.point1.y - edge2.point2.y;
        PointType b2 = edge2.point2.x - edge2.point1.x;
        PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

        PointType line1_start = a2 * edge1.point1.x + b2 * edge1.point1.y + d2;
        PointType line1_end   = a2 * edge1.point2.x + b2 * edge1.point2.y + d2;

        PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
        PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

        PointType u         = line1_start / (line1_start - line1_end);
        out_intersection->x = start1.x + u * dir11.x;
        out_intersection->y = start1.y + u * dir11.y;
        out_intersection->z = start1.z + u * dir11.z;

        if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
            return true;

        return false;
    }

    template <class PointType> int IsOnEdge(DGeo::Point<PointType> point, PointType epsilon)
    {
        if (std::abs((point.x - point1.x) * (point2.y - point1.y) - (point.y - point1.y) * (point2.x - point1.x)) < epsilon)
        {
            if ((point1.x <= point.x && point2.x > point.x) || (point2.x <= point.x && point1.x > point.x) ||
                (point1.y <= point.y && point2.y > point.y) || (point2.y <= point.y && point1.y > point.y))
            {
                if (point1.y == point2.y)
                    return 1;
                return 0;
            }
        };
        return 2;
    }
};
}
#endif