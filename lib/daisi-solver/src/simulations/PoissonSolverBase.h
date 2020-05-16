#ifndef POISSONSOLVERBASE_H
#define POISSONSOLVERBASE_H
#include "Geom.h"
namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
};

template <class PointType>
class BoundaryPoint
{
  public:
    int                   pnum;  // number of mesh point
    char                  type;  // type of boundary condition
    double                value; // value of boundary condition
    DGeo::Edge<PointType> edge;
    int                   bcnum; // number of boundary condition;
};

template <class PointType>
class IBPoint
{ // connnection between intersection point and mesh boundary point
  public:
    int                    bpnum;  // number of boundary point
    DGeo::Point<PointType> ipoint; // coordinates of intersection point
};

template <class PointType>
class BoundaryPointIntersection
{
  public:
    char                   type;
    double                 value;
    DGeo::Edge<PointType>  edge;
    DGeo::Point<PointType> point;
    int                    bcnum;
};

template <class PointType>
class NonZeroLinVect
{ // non zero points in linVect;
  public:
    int       row;
    PointType coef;
    int       z;
    int       bcnum;
};
#endif