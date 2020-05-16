#include "Geom.h"
#include "VTKIncludeSolver.h"
#include <common_tools/constants.h>

using namespace DGeo;

template class Point<float>;
template class Point<double>;
template class Point<int>;
template class Edge<float>;
template class Edge<double>;

template <class PointType>
bool Point<PointType>::isEqual(Point<PointType> point1, PointType epsilon) const
{
    if (std::abs(point1.x - x) < epsilon && std::abs(point1.y - y) < epsilon && std::abs(point1.z - z) < epsilon)
        return true;
    return false;
};

template <class PointType>
PointType Point<PointType>::Radius()
{
    return sqrt(x * x + y * y + z * z);
};
template <class PointType>
PointType Point<PointType>::Dist2Point(const Point<PointType>& point1) const
{
    return sqrt((x - point1.x) * (x - point1.x) + (y - point1.y) * (y - point1.y) +
                (z - point1.z) * (z - point1.z));
};
template <class PointType>
Point<PointType> Point<PointType>::rotate(PointType alpha, Point<PointType> point)
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
template <class PointType>
Point<PointType>& Point<PointType>::operator=(const Point<PointType>& right)
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

template <class PointType>
Point<PointType> Point<PointType>::operator+(const Point<PointType>& right) const
{
    Point<PointType> result = *this;
    result.x += right.x;
    result.y += right.y;
    result.z += right.z;
    return result;
};
template <class PointType>
Point<PointType> Point<PointType>::operator-(const Point<PointType>& right) const
{
    Point<PointType> result = *this;
    result.x -= right.x;
    result.y -= right.y;
    result.z -= right.z;
    return result;
};
template <class PointType>
Point<PointType>::Point()
{
    isOut = false;
};
template <class PointType>
Point<PointType>::Point(Point<PointType> Old, PointType hx, PointType hy, PointType hz)
{
    *this = Old;
    x += hx;
    y += hy;
    z += hz;
};

template <class PointType>
bool Edge<PointType>::IsOnOneSide(const DGeo::Point<PointType>& pp1,
                                  const DGeo::Point<PointType>& pp2, double epsilon) const
{
    epsilon = epsilon * length();

    PointType a2 = point1.y - point2.y;
    PointType b2 = point2.x - point1.x;
    PointType d2 = point1.x * point2.y - point2.x * point1.y;

    PointType line1_start = a2 * pp1.x + b2 * pp1.y + d2;
    PointType line1_end   = a2 * pp2.x + b2 * pp2.y + d2;
    if (line1_start * line1_end < -epsilon * epsilon)
        return false;
    return true;
}

template <class PointType>
Point<PointType> Edge<PointType>::GetNormalPoint1(PointType h) const
{
    double           al  = alpha();
    Point<PointType> mid = Middle();
    Point<PointType> result;

    result.x = mid.x + h * std::sin(al);
    result.y = mid.y - h * std::cos(al);
    result.z = mid.z;

    return result;
}
template <class PointType>
Point<PointType> Edge<PointType>::GetNormalPoint2(PointType h) const
{
    double           al  = alpha();
    Point<PointType> mid = Middle();
    Point<PointType> result;

    result.x = mid.x - h * std::sin(al);
    result.y = mid.y + h * std::cos(al);
    result.z = mid.z;

    return result;
}

template <class PointType>
vtkSmartPointer<vtkUnstructuredGrid> Edge<PointType>::GetVTKGrid()
{
    vtkSmartPointer<vtkUnstructuredGrid> VTKgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>           points  = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(point1.x, point1.y, point1.z);
    points->InsertNextPoint(point2.x, point2.y, point2.z);

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 0);
    line->GetPointIds()->SetId(1, 1);
    cellArray->InsertNextCell(line);

    VTKgrid->SetPoints(points);
    VTKgrid->SetCells(VTK_LINE, cellArray);

    return VTKgrid;
};

template <class PointType>
PointType Edge<PointType>::alpha() const
{
    if (point1.x == point2.x)
        return commtools::PI() / 2;

    /*	r = sqrt(x*x + y*y);
            if (x > 0 && y >= 0)
            {
                    return phi = atan(y / x);
            }
            if (x > 0 && y < 0)
            {
                    phi = atan(y / x) + 2 * commtools::PI();

            }
            if (x < 0)
            {
                    phi = atan(y / x) + commtools::PI();

            }
            if (std::abs(x)<1e-7 && y>0)
            {
                    phi = commtools::PI() / 2;

            }
            if (std::abs(x)<1e-7 && y<0)
            {
                    phi = 3 * commtools::PI() / 2;

            }*/

    // return atan(std::abs((point1.y - point2.y) / ((point1.x - point2.x))));

    PointType y = point1.y - point2.y;
    PointType x = point1.x - point2.x;
    PointType l = length();

    if (y >= 0)
        return std::acos(x / l);
    else
        return -std::acos(x / l);

    // return atan(((point1.y - point2.y) / ((point1.x - point2.x))));
}
template <class PointType>
Point<PointType> Edge<PointType>::Middle() const
{
    Point<PointType> p;
    p.x = (point1.x + point2.x) / 2;
    p.y = (point1.y + point2.y) / 2;
    p.z = (point1.z + point2.z) / 2;
    return p;
};
template <class PointType>
std::vector<Edge<PointType>> Edge<PointType>::resize(int size)
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
template <class PointType>
Edge<PointType>& Edge<PointType>::operator=(const Edge<PointType>& right)
{
    if (this == &right)
    {
        return *this;
    }
    point1 = right.point1;
    point2 = right.point2;
    return *this;
};
template <class PointType>
PointType Edge<PointType>::length() const
{
    return sqrt((point1.x - point2.x) * (point1.x - point2.x) +
                (point1.y - point2.y) * (point1.y - point2.y) +
                (point1.z - point2.z) * (point1.z - point2.z));
};
template <class PointType>
bool Edge<PointType>::IsParallel(const DGeo::Edge<PointType>& edge2, double epsilon) const
{
    //	epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

    epsilon = std::min(epsilon * length(), edge2.length() * epsilon);

    /*float k1 = std::abs((edge2.point2.y - edge2.point1.y) / ((edge2.point2.x - edge2.point1.x)));

    float k2 = std::abs((point2.y - edge2.point1.y) / ((point2.x - point1.x)));

    if (std::abs(k1 - k1) - epsilon)
            return true;

    if ((edge2.point2.x - edge2.point1.x)<epsilon && (point2.x - point1.x)<epsilon)
            return true;*/

    /*if ((edge2.point2.x - edge2.point1.x)<epsilon && (point2.x - point1.x)>epsilon)
            return false;

    if ((edge2.point2.x - edge2.point1.x)>epsilon && (point2.x - point1.x)<epsilon)
            return false;

    if ((edge2.point2.y - edge2.point1.y)<epsilon && (point2.y - point1.y)>epsilon)
            return false;

    if ((edge2.point2.y - edge2.point1.y)>epsilon && (point2.y - point1.y)<epsilon)
            return false;*/

    if (std::abs((edge2.point1.x - point1.x) * (point2.y - point1.y) -
            (edge2.point1.y - point1.y) * (point2.x - point1.x)) < epsilon &&
        std::abs((edge2.point2.x - point1.x) * (point2.y - point1.y) -
            (edge2.point2.y - point1.y) * (point2.x - point1.x)) < epsilon)
        return true;

    return false;
};
template <class PointType>
bool Edge<PointType>::IsParallelFixEps(const DGeo::Edge<PointType>& edge2, double epsilon) const
{
    /*if ((edge2.point2.x - edge2.point1.x)<epsilon && (point2.x - point1.x)>epsilon)
            return false;

    if ((edge2.point2.x - edge2.point1.x)>epsilon && (point2.x - point1.x)<epsilon)
            return false;

    if ((edge2.point2.y - edge2.point1.y)<epsilon && (point2.y - point1.y)>epsilon)
            return false;

    if ((edge2.point2.y - edge2.point1.y)>epsilon && (point2.y - point1.y)<epsilon)
            return false;*/

    if (std::abs((edge2.point1.x - point1.x) * (point2.y - point1.y) -
            (edge2.point1.y - point1.y) * (point2.x - point1.x)) < epsilon &&
        std::abs((edge2.point2.x - point1.x) * (point2.y - point1.y) -
            (edge2.point2.y - point1.y) * (point2.x - point1.x)) < epsilon)
        return true;

    return false;
};
template <class PointType>
bool Edge<PointType>::IsIntersectionEdge(const DGeo::Edge<PointType>& edge2, double epsilon) const
{
    // epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

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

    epsilon = std::min(epsilon * length(), edge2.length() * epsilon);

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

template <class PointType>
bool Edge<PointType>::IsIntersectionEdgeMesh(const DGeo::Edge<PointType>& edge2,
                                             double                       epsilon) const
{
    // epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

    /*DGeo::Edge<PointType> edge1 = *this;
    DGeo::Point<PointType> start1 = edge1.point1;
    DGeo::Point<PointType> start2 = edge2.point1;
    DGeo::Point<PointType> end1 = edge1.point2;
    DGeo::Point<PointType> end2 = edge2.point2;

    if (edge1.IsParallel(edge2, epsilon))
            return false;*/

    // epsilon = std::min(epsilon*length(), edge2.length()*epsilon);

    //������� ��������� ������ ���������� ����� �������
    PointType a1 = point1.y - point2.y;
    PointType b1 = point2.x - point1.x;
    PointType d1 = point1.x * point2.y - point2.x * point1.y;

    PointType a2 = edge2.point1.y - edge2.point2.y;
    PointType b2 = edge2.point2.x - edge2.point1.x;
    PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

    //����������� ����� ��������, ��� ��������� � ����� ������������� ���
    PointType line1_start = a2 * point1.x + b2 * point1.y + d2;
    PointType line1_end   = a2 * point2.x + b2 * point2.y + d2;

    PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
    PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

    //���� ����� ������ ������� ����� ���� ����, ������ �� � ����� ������������� � ����������� ���.
    if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
        return true;

    //	PointType u = seg1_line2_start / (seg1_line2_start - seg1_line2_end);
    //	*out_intersection = start1 + u*dir1;

    return false;
}

template <class PointType>
bool Edge<PointType>::IsIntersectionEdge(const DGeo::Edge<PointType>& edge2, double epsilon,
                                         DGeo::Point<PointType>* out_intersection) const
{
    // epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

    //

    if (IsParallel(edge2, epsilon))
        return false;

    if (IsOnEdge(edge2.point2, epsilon) < 2)
    {
        out_intersection->x = edge2.point2.x;
        out_intersection->y = edge2.point2.y;
        out_intersection->z = 0;
        return true;
    }
    if (IsOnEdge(edge2.point1, epsilon) < 2)
    {
        out_intersection->x = edge2.point1.x;
        out_intersection->y = edge2.point1.y;
        out_intersection->z = 0;

        return true;
    }

    if (edge2.IsOnEdge(point2, epsilon) < 2)
    {
        out_intersection->x = point2.x;
        out_intersection->y = point2.y;
        out_intersection->z = 0;

        return true;
    }
    if (edge2.IsOnEdge(point1, epsilon) < 2)
    {
        out_intersection->x = point1.x;
        out_intersection->y = point1.y;
        out_intersection->z = 0;

        return true;
    }

    epsilon = std::min(epsilon * length(), edge2.length() * epsilon);

    PointType a1 = point1.y - point2.y;
    PointType b1 = point2.x - point1.x;
    PointType d1 = point1.x * point2.y - point2.x * point1.y;

    PointType a2 = edge2.point1.y - edge2.point2.y;
    PointType b2 = edge2.point2.x - edge2.point1.x;
    PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

    PointType line1_start = a2 * point1.x + b2 * point1.y + d2;
    PointType line1_end   = a2 * point2.x + b2 * point2.y + d2;

    PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
    PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

    PointType u         = line1_start / (line1_start - line1_end);
    out_intersection->x = point1.x + u * b1;
    out_intersection->y = point1.y - u * a1;
    out_intersection->z = point1.z + u * (point2.z - point1.z);

    if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
        return true;

    return false;
}
template <class PointType>
bool Edge<PointType>::IsIntersectionEdgeLight(const DGeo::Edge<PointType>& edge2, double epsilon,
                                              DGeo::Point<PointType>* out_intersection) const
{
    //	epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));
    // epsilon = std::min(epsilon*length(), edge2.length()*epsilon);

    PointType a1 = point1.y - point2.y;
    PointType b1 = point2.x - point1.x;
    PointType d1 = point1.x * point2.y - point2.x * point1.y;

    PointType a2 = edge2.point1.y - edge2.point2.y;
    PointType b2 = edge2.point2.x - edge2.point1.x;
    PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

    PointType line1_start = a2 * point1.x + b2 * point1.y + d2;
    PointType line1_end   = a2 * point2.x + b2 * point2.y + d2;

    PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
    PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

    PointType u         = line1_start / (line1_start - line1_end);
    out_intersection->x = point1.x + u * b1;
    out_intersection->y = point1.y - u * a1;
    out_intersection->z = point1.z + u * (point2.z - point1.z);

    // if (line1_start * line1_end < -pow(double(epsilon), 1.5) && line2_start * line2_end <
    // -pow(double(epsilon), 1.5))

    //	if ((line1_start * line1_end < 0 || std::abs(line1_start * line1_end)<epsilon*epsilon) &&
    //(line2_start * line2_end <
    // 0 || std::abs(line2_start * line2_end)<epsilon*epsilon)) 		return true;

    if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
        return true;

    return false;
}
template <class PointType>
bool Edge<PointType>::IsIntersectionEdgePoints(const DGeo::Edge<PointType>& edge2, double epsilon,
                                               DGeo::Point<PointType>* out_intersection) const
{
    // epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

    if (IsParallel(edge2, epsilon))
        return false;

    if (IsOnEdge(edge2.point2, epsilon) < 2)
    {
        out_intersection->x = edge2.point2.x;
        out_intersection->y = edge2.point2.y;
        return true;
    }
    if (IsOnEdge(edge2.point1, epsilon) < 2)
    {
        out_intersection->x = edge2.point1.x;
        out_intersection->y = edge2.point1.y;
        return true;
    }

    if (edge2.IsOnEdge(point2, epsilon) < 2)
    {
        out_intersection->x = point2.x;
        out_intersection->y = point2.y;
        return true;
    }
    if (edge2.IsOnEdge(point1, epsilon) < 2)
    {
        out_intersection->x = point1.x;
        out_intersection->y = point1.y;
        return true;
    }

    epsilon = std::min(epsilon * length(), edge2.length() * epsilon);

    PointType a1 = point1.y - point2.y;
    PointType b1 = point2.x - point1.x;
    PointType d1 = point1.x * point2.y - point2.x * point1.y;

    PointType a2 = edge2.point1.y - edge2.point2.y;
    PointType b2 = edge2.point2.x - edge2.point1.x;
    PointType d2 = edge2.point1.x * edge2.point2.y - edge2.point2.x * edge2.point1.y;

    PointType line1_start = a2 * point1.x + b2 * point1.y + d2;
    PointType line1_end   = a2 * point2.x + b2 * point2.y + d2;

    PointType line2_start = a1 * edge2.point1.x + b1 * edge2.point1.y + d1;
    PointType line2_end   = a1 * edge2.point2.x + b1 * edge2.point2.y + d1;

    PointType u         = line1_start / (line1_start - line1_end);
    out_intersection->x = point1.x + u * b1;
    out_intersection->y = point1.y - u * a1;
    out_intersection->z = point1.z + u * (point2.z - point1.z);

    if (line1_start * line1_end < 0 && line2_start * line2_end < 0)
        return true;

    return false;
}

template <class PointType>
int Edge<PointType>::IsOnEdge(const DGeo::Point<PointType>& point, double epsilon) const
{
    //	epsilon = std::min(double(Dconst::epsFrac*length()), double(epsilon));

    // if (point.isEqual(point1, epsilon) || point.isEqual(point2, epsilon))
    //	return 1;

    epsilon = epsilon * length();
    double E = 1e-16;
    if (std::abs((point.x - point1.x) * (point2.y - point1.y) -
            (point.y - point1.y) * (point2.x - point1.x)) < epsilon)
    {
        PointType l1 = sqrt((point1.x - point.x) * (point1.x - point.x) +
                            (point1.y - point.y) * (point1.y - point.y));
        PointType l2 = sqrt((point2.x - point.x) * (point2.x - point.x) +
                            (point2.y - point.y) * (point2.y - point.y));
        PointType l = length();
        //	if (std::abs(std::abs(point1.x - point.x) + std::abs(point2.x - point.x) - std::abs(point1.x -
        // point2.x))<200 * epsilon &&
        // std::abs(std::abs(point1.y - point.y) + std::abs(point2.y - point.y) - std::abs(point1.y - point2.y))<200 *
        // epsilon)

//        if ((((point.x < point1.x && point.x > point2.x) ||
//              (point.x > point1.x && point.x < point2.x) || std::abs(point.x - point1.x)<E ||
//                std::abs(point.x - point2.x)<E) &&
//             ((point.y < point1.y && point.y > point2.y) ||
//              (point.y > point1.y && point.y < point2.y) || std::abs (point.y - point1.y)<E ||
//              std::abs(point.y - point2.y)<E)) ||
//            std::abs(l1 + l2 - l) < epsilon)

          if ((((point.x < point1.x && point.x > point2.x) ||
                (point.x > point1.x && point.x < point2.x) || point.x == point1.x ||
                point.x == point2.x) &&
               ((point.y < point1.y && point.y > point2.y) ||
                (point.y > point1.y && point.y < point2.y) || point.y == point1.y ||
                point.y == point2.y)) ||
              std::abs(l1 + l2 - l) < epsilon)
        //(std::abs(point.x - point1.x)<epsilon ) || std::abs(point.x - point2.x)<epsilon || std::abs(point.y -
        // point1.y)<epsilon ||
        // std::abs(point.y - point2.y)<epsilon)
        {
            //if ( std::abs(point1.y - point2.y) <E)
            if ( point1.y == point2.y)
                return 1;
            return 0;
        }
        /*if (std::abs(l1 + l2 - l)< 10 * epsilon)
        {
                if (point1.y == point2.y)
                        return 1;
                return 0;
        }*/
    };
    return 2;
}

template <class PointType>
void Edge<PointType>::SwapPoints()
{
    Point<PointType> tmp = point1;
    point1               = point2;
    point2               = tmp;
};

template <class PointType>
bool Edge<PointType>::IsOnLine(const DGeo::Point<PointType>& point, double epsilon) const
{
    // epsilon = 2 * std::min(double(Dconst::epsFrac*length()), double(epsilon));

    epsilon = epsilon * length();

    if (point.isEqual(point1, epsilon) || point.isEqual(point2, epsilon))
        return true;

    if (std::abs((point.x - point1.x) * (point2.y - point1.y) -
            (point.y - point1.y) * (point2.x - point1.x)) < epsilon)
    {
        return true;
    };
    return false;
}

template <class PointType>
bool Edge<PointType>::IsEqual(const DGeo::Edge<PointType>& edge2, double epsilon) const
{
    epsilon = std::min(epsilon * length(), edge2.length() * epsilon);
    if ((point1.isEqual(edge2.point1, epsilon) && point2.isEqual(edge2.point2, epsilon)) ||
        (point1.isEqual(edge2.point2, epsilon) && point2.isEqual(edge2.point1, epsilon)))
        return true;
    return false;
};

template <class PointType>
Edge<PointType> Edge<PointType>::TranslateByNormal1(PointType h) const
{
    double          al     = alpha();
    Edge<PointType> result = *this;

    result.point1.x = result.point1.x + h * std::sin(al);
    result.point1.y = result.point1.y - h * std::cos(al);
    result.point1.z = result.point1.z;

    result.point2.x = result.point2.x + h * std::sin(al);
    result.point2.y = result.point2.y - h * std::cos(al);
    result.point2.z = result.point2.z;

    return result;
};

template <class PointType>
Edge<PointType> Edge<PointType>::TranslateByNormal2(PointType h) const
{
    double          al     = alpha();
    Edge<PointType> result = *this;

    result.point1.x = result.point1.x - h * std::sin(al);
    result.point1.y = result.point1.y + h * std::cos(al);
    result.point1.z = result.point1.z;

    result.point2.x = result.point2.x - h * std::sin(al);
    result.point2.y = result.point2.y + h * std::cos(al);
    result.point2.z = result.point2.z;

    return result;
};