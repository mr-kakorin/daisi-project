#include "PoissonSolver3d.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "GridData.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"

//#include "SOR.h"

//#include "mkl_vml.h"

template class PoissonSolver3d<float>;
template class PoissonSolver3d<double>;

template <class PointType>
void PoissonSolver3d<PointType>::SetParameters(const std::vector<double>& param)
{
    w                    = param[0];
    eps_tolerance        = param[1];
    w_charge             = param[2];
    eps_tolerance_charge = param[3];

    eps_p = param[4];
    eps_h = param[5];

    ChargeSpace          = param[6];
    RecalculateParameter = param[7];

    for (int i         = 0; i < solverFlags.size(); i++)
        solverFlags[i] = param[8 + i];
};

template <class PointType>
std::vector<double> PoissonSolver3d<PointType>::GetParameters()
{
    std::vector<double> result(8);

    result[0] = w;
    result[1] = eps_tolerance;

    result[2] = w_charge;
    result[3] = eps_tolerance_charge;

    result[4] = eps_p;
    result[5] = eps_h;

    result[6] = ChargeSpace;
    result[7] = RecalculateParameter;

    for (int i = 0; i < solverFlags.size(); i++)
        result.push_back(int(solverFlags[i]));

    return result;
};

template <class PointType>
bool compNum(BoundaryPoint3d<PointType> a, BoundaryPoint3d<PointType> b)
{
    return (a.pnum < b.pnum);
};

int sign(double val)
{
    if (val == 0)
        return 0;
    if (val < 0)
        return -1;
    return 1;
}

// hx1 - left, hx2 - right, hy1 - down, hy2 - up;
double deep3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
              double y, double z)
{
    return 2 / (hz2 * (hz1 + hz2)); // poisson
}

double outward3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
                 double y, double z)
{
    return 2 / (hz1 * (hz2 + hz1)); // poisson
}

double up3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
            double y, double z)
{
    return 2 / (hy2 * (hy1 + hy2)); // poisson
}

double down3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
              double y, double z)
{
    return 2 / (hy1 * (hy2 + hy1)); // poisson
}

double left3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
              double y, double z)
{
    return 2 / (hx1 * (hx1 + hx2)); // poisson
}

double right3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
               double y, double z)
{
    return 2 / (hx2 * (hx1 + hx2)); // poisson
}

double middle3d(double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
                double y, double z)
{
    return -2 / (hx1 * hx2) - 2 / (hy1 * hy2) - 2 / (hz1 * hz2); // poisson
}

template <class PointType>
void PoissonSolver3d<PointType>::searchBoundaryPointsIntersection_3d(
    const std::shared_ptr<MeshContainer3d<PointType>>&           mesh,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaries,
    const std::shared_ptr<BoundaryConditions>& boundaryConditions, double z)
{
    boundaryPointsIntersectionCondition.clear();
    boundaryPointsIntersectionDefault.clear();
    DGeo::Edge<PointType> tmpEdge;
    tmpEdge.point1.z = 0;
    tmpEdge.point2.z = 0;
    std::vector<int>                    boundariesList;
    std::vector<DGeo::Point<PointType>> intersectionPointsTmp;
    std::vector<DGeo::Edge<PointType>>  intersectionEdgesTmp;

    BoundaryPointIntersection<PointType> pointTmp;
    boundCondNum = boundaryConditions->PropertyConditionListSize();

    int status;

    for (int i = 0; i < mesh->mesh[0].meshData.size(); i++)
    {

        tmpEdge.point1.x = mesh->mesh[0].meshData[i][0].x;
        tmpEdge.point2.x = mesh->mesh[0].meshData[i].back().x;

        tmpEdge.point1.y = mesh->mesh[0].meshData[i][0].y;
        tmpEdge.point2.y = mesh->mesh[0].meshData[i].back().y;

        for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++)
        {
            boundaryConditions->GetPotential(j, 0, z, status);
            if (status == 1)
            {
                boundariesList = boundaryConditions->GetPropertyConditionsBoundariesList(j);
                for (int k = 0; k < boundariesList.size(); k++)
                {
                    boundaries[boundariesList[k]]->FindAllIntersections(
                        eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp);
                    for (int s = 0; s < intersectionPointsTmp.size(); s++)
                    {
                        pointTmp.point = intersectionPointsTmp[s];
                        pointTmp.edge  = intersectionEdgesTmp[s];
                        pointTmp.type  = 'd';
                        pointTmp.value = boundaryConditions->GetPotentialOffset(j);
                        pointTmp.bcnum = j;
                        boundaryPointsIntersectionCondition.push_back(pointTmp);
                    };
                }
            }
        };
    }

    for (int i = 0; i < mesh->mesh[0].getnVertexX() - 1; i++)
    {
        for (int j = 0; j < mesh->mesh[0].getnVertexY() - 1; j++)
        {
            if (mesh->templNumb(i, j, 1) != -1)
            {
                if (mesh->templNumb(i, j + 1, 1) != -1)
                {
                    if (mesh->serialMeshData[mesh->templNumb(i, j, 1)].isOut == false &&
                        mesh->serialMeshData[mesh->templNumb(i, j + 1, 1)].isOut == false)
                        continue;

                    tmpEdge.point1.x = mesh->serialMeshData[mesh->templNumb(i, j, 1)].x;
                    tmpEdge.point2.x = mesh->serialMeshData[mesh->templNumb(i, j + 1, 1)].x;

                    tmpEdge.point1.y = mesh->serialMeshData[mesh->templNumb(i, j, 1)].y;
                    tmpEdge.point2.y = mesh->serialMeshData[mesh->templNumb(i, j + 1, 1)].y;

                    for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++)
                    {
                        boundaryConditions->GetPotential(j, 0, z, status);
                        if (status == 1)
                        {
                            boundariesList =
                                boundaryConditions->GetPropertyConditionsBoundariesList(j);
                            for (int k = 0; k < boundariesList.size(); k++)
                            {
                                boundaries[boundariesList[k]]->FindAllIntersections(
                                    eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp);

                                for (int s = 0; s < intersectionPointsTmp.size(); s++)
                                {
                                    pointTmp.point = intersectionPointsTmp[s];
                                    pointTmp.edge  = intersectionEdgesTmp[s];
                                    pointTmp.type  = 'd';
                                    pointTmp.value = boundaryConditions->GetPotentialOffset(j);
                                    pointTmp.bcnum = j;
                                    boundaryPointsIntersectionCondition.push_back(pointTmp);
                                };
                            }
                        }
                    }
                }
            }
        }
    }
};

template <class PointType>
DGeo::Point<int> PoissonSolver3d<PointType>::closestVertex_3d(
    const std::shared_ptr<MeshContainer3d<PointType>>& mesh, DGeo::Point<PointType> point)
{
    DGeo::Point<int> result;
    int              i, j, k;
    i        = 1; // x
    j        = 1; // y
    k        = 1; // z
    int flag = 0;
    // z
    while (mesh->templNumb(i, j, k) == -1)
    {
        if (j + 1 > mesh->templNumb.GetNcol() - 1)
        {
            ++i;
            j = 1;
        }
        else
        {
            ++j;
        }
    }

    while ((mesh->serialMeshData[mesh->templNumb(i, j, k)].z <= point.z) &&
           (k + 1 < mesh->templNumb.GetNz() - 1))
    {
        ++k;
        i = 1;
        j = 1;
        while (mesh->templNumb(i, j, k) == -1)
        {
            if (j + 1 > mesh->templNumb.GetNcol() - 1)
            {
                ++i;
                j = 1;
            }
            else
            {
                ++j;
            }
        }
    }
    if (std::abs(point.z - mesh->serialMeshData[mesh->templNumb(i, j, k)].z) >
        std::abs(point.z - mesh->serialMeshData[mesh->templNumb(i, j, k - 1)].z))
    {
        --k;
    }
    // y

    i = 1;
    j = 1;
    while (mesh->templNumb(i, j, k) == -1)
    {
        if (i + 1 < mesh->templNumb.GetNrow() - 1)
        {
            ++i;
        }
        else
        {
            ++j;
            i = 1;
        }
    }

    while ((mesh->serialMeshData[mesh->templNumb(i, j, k)].y < point.y) &&
           (j + 1 < mesh->templNumb.GetNcol() - 1))
    {
        double y = mesh->serialMeshData[mesh->templNumb(i, j, k)].y;
        i        = 1;
        ++j;
        while (mesh->templNumb(i, j, k) == -1)
        {
            if (i + 1 < mesh->templNumb.GetNrow() - 1)
            {
                ++i;
            }
            else
            {
                ++j;
                i = 1;
            }
        }
    }

    if (std::abs(point.y - mesh->serialMeshData[mesh->templNumb(i, j, k)].y) >
        std::abs(point.y - mesh->serialMeshData[mesh->templNumb(i, j - 1, k)].y))
    {
        --j;
    }

    // x
    i = 1;
    while (mesh->templNumb(i, j, k) == -1)
    {
        ++i;
    }
    while ((mesh->serialMeshData[mesh->templNumb(i, j, k)].x < point.x) &&
           (i + 1 < mesh->templNumb.GetNrow() - 1))
    {
        ++i;
    }
    if (std::abs(point.y - mesh->serialMeshData[mesh->templNumb(i, j, k)].y) >
        std::abs(point.y - mesh->serialMeshData[mesh->templNumb(i - 1, j, k)].y))
    {
        --i;
    }

    result.x = i;
    result.y = j;
    result.z = k;

    return result;
};

template <class PointType>
double PoissonSolver3d<PointType>::dist(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2)
{

    double dx, dy;
    dx = vert2.x - vert1.x;
    dy = vert2.y - vert1.y;

    double result = sqrt(dx * dx + dy * dy);
    return result;
};

template <class PointType>
double PoissonSolver3d<PointType>::distToEdgeX(DGeo::Edge<PointType>  edge,
                                               DGeo::Point<PointType> point)
{

    if (polar == false)
    {
        double y1, y2, x1, x2;
        y1 = edge.point1.y;
        y2 = edge.point2.y;
        x1 = edge.point1.x;
        x2 = edge.point2.x;

        if (edge.IsOnLine(point, eps_p))
        {
            double h1 = sqrt((x1 - point.x) * (x1 - point.x) + (y1 - point.y) * (y1 - point.y));
            double h2 = sqrt((x2 - point.x) * (x2 - point.x) + (y2 - point.y) * (y2 - point.y));

            return PointType(std::min(h1, h2));
        }
        if (std::abs(x1 - x2) < eps_p)
            return std::abs(x1 - point.x);
        else
        {
            double k = (y2 - y1) / (x2 - x1);
            double b = y1 - k * x1;
            return std::abs(point.x - (point.y - b) / k);
        }
        return 0;
    }
    else
    {
        PointType pi = 3.14159265358979;
        PointType k1, k2;
        PointType b;

        PointType y1, y2, x1, x2;
        y1 = edge.point1.y;
        y2 = edge.point2.y;
        x1 = edge.point1.x;
        x2 = edge.point2.x;

        PointType x, y, dx, dy;
        Dmath::Polar2Cartesian(point.x, point.y, x, y);

        if (std::abs(y1 - y2) < eps_p)
        {
            dx = 0;
            dy = y - y1;
        }
        else if (std::abs(x1 - x2) < eps_p)
        {
            dx = x - x1;
            dy = 0;
        }
        else
        {
            k1 = (y2 - y1) / (x2 - x1);
            b  = y1 - k1 * x1;

            k2 = y / x;
            dx = x - b / (k2 - k1);
            dy = y - k2 * b / (k2 - k1);
        }
        return std::abs(sqrt(dx * dx + dy * dy));
    }
    return 0;
};

template <class PointType>
double PoissonSolver3d<PointType>::distToEdgeY(DGeo::Edge<PointType>  edge,
                                               DGeo::Point<PointType> point)
{

    if (polar == false)
    {
        double y1, y2, x1, x2;
        y1 = edge.point1.y;
        y2 = edge.point2.y;
        x1 = edge.point1.x;
        x2 = edge.point2.x;

        if (edge.IsOnLine(point, eps_p))
        {
            double h1 = sqrt((x1 - point.x) * (x1 - point.x) + (y1 - point.y) * (y1 - point.y));
            double h2 = sqrt((x2 - point.x) * (x2 - point.x) + (y2 - point.y) * (y2 - point.y));

            return PointType(std::min(h1, h2));
        }

        if (std::abs(y1 - y2) < eps_p)
            return std::abs(y1 - point.y);
        else
        {
            double k = (y2 - y1) / (x2 - x1);
            double b = y1 - k * x1;
            return std::abs(point.y - k * point.x - b);
        }
    }
    else
    {
        PointType pi = 3.14159265358979;
        PointType k;
        PointType b;

        PointType y1, y2, x1, x2;
        y1 = edge.point1.y;
        y2 = edge.point2.y;
        x1 = edge.point1.x;
        x2 = edge.point2.x;

        PointType x, y, r;
        PointType x_1, x_2, y_1, y_2;
        double    dphi_1, dphi_2;
        Dmath::Polar2Cartesian(point.x, point.y, x, y);
        r = sqrt(x * x + y * y);
        if (std::abs(y1 - y2) < eps_p)
        {
            if (r * r < y1 * y1)
            {
                return 0;
            }
            x_1 = sqrt(r * r - y1 * y1);
            x_2 = -sqrt(r * r - y1 * y1);
            y_1 = y1;
            y_2 = y1;
        }
        else if (std::abs(x1 - x2) < eps_p)
        {
            if (r * r < x1 * x1)
            {
                return 0;
            }
            x_1 = x1;
            x_2 = x1;
            y_1 = sqrt(r * r - x1 * x1);
            y_2 = -sqrt(r * r - x1 * x1);
        }
        else
        {
            k = (y2 - y1) / (x2 - x1);
            b = y1 - k * x1;

            x_1 = (-2 * k * b + sqrt(4 * k * k * b * b - (4 * (b * b - r * r)) * (k + 1))) /
                  (2 * (k + 1));
            x_2 = (-2 * k * b - sqrt(4 * k * k * b * b - (4 * (b * b - r * r)) * (k + 1))) /
                  (2 * (k + 1));

            y_1 = k * x_1 + b;
            y_2 = k * x_2 + b;
        }
        PointType r_1, r_2, phi_1, phi_2;

        Dmath::Cartesian2Polar(x_1, y_1, r_1, phi_1);
        Dmath::Cartesian2Polar(x_2, y_2, r_2, phi_2);

        dphi_1 = std::abs(phi_1 - point.y);
        dphi_2 = std::abs(phi_2 - point.y);

        if (dphi_1 < dphi_2)
        {
            return dphi_1;
        }
        else
        {
            return dphi_2;
        }
    }
    return 0;
};

template <class PointType>
double PoissonSolver3d<PointType>::distP2PX(DGeo::Point<PointType> vert1,
                                            DGeo::Point<PointType> vert2)
{
    return std::abs(vert2.x - vert1.x);
};

template <class PointType>
double PoissonSolver3d<PointType>::distP2PY(DGeo::Point<PointType> vert1,
                                            DGeo::Point<PointType> vert2)
{
    return std::abs(vert2.y - vert1.y);
};

template <class PointType>
double PoissonSolver3d<PointType>::distP2PZ(DGeo::Point<PointType> vert1,
                                            DGeo::Point<PointType> vert2)
{
    return std::abs(vert2.z - vert1.z);
};

template <class PointType>
double PoissonSolver3d<PointType>::distP2P(DGeo::Point<PointType> vert1,
                                           DGeo::Point<PointType> vert2, char type)
{
    if (type == 'x')
    {
        return std::abs(vert2.x - vert1.x);
    }
    else if (type == 'y')
    {
        return std::abs(vert2.y - vert1.y);
    }
    else
    {
        return std::abs(vert2.z - vert1.z);
    }
};

template <class PointType>
DGeo::Point<PointType> normalVectToEdge(DGeo::Edge<PointType> edge)
{

    DGeo::Point<PointType> result;
    double                 x0, y0, x, y;

    x0 = edge.point2.x - edge.point1.x;
    y0 = edge.point2.y - edge.point2.y;

    x = 1;
    y = -x0 / y0;

    double normk = 1 / (sqrt(x * x + y * y));
    result.x     = x * normk;
    result.y     = y * normk;
    return result;
};

template <class PointType>
BoundaryPoint3d<PointType> PoissonSolver3d<PointType>::findBoundaryPoint(DGeo::Point<int> bpoint)
{

    int i = 0;
    if (boundaryPoints[bpoint.z - 1][bpoint.y - 1].size() > 0)
    {
        for (i = 0; i < boundaryPoints[bpoint.z - 1][bpoint.y - 1].size(); ++i)
        {
            if (bpoint.Number == boundaryPoints[bpoint.z - 1][bpoint.y - 1][i].pnum)
            {
                return boundaryPoints[bpoint.z - 1][bpoint.y - 1][i];
            }
        }
    }
    BoundaryPoint3d<PointType> bp;
    bp.type  = 'f';
    bp.value = 0;
    bp.pnum  = -1;
    bp.bcnum = -2;

    return bp;
};

template <class PointType>
void PoissonSolver3d<PointType>::addNearBoundaryPoint_n(int n)
{

    if (nearBoundaryPoints_n.size() == 0)
    {
        nearBoundaryPoints_n.push_back(n);
        return;
    }
    int i = 0;
    for (i = 0; i < nearBoundaryPoints_n.size(); ++i)
    {
        if (n >= nearBoundaryPoints_n[i])
        {
            break;
        }
    }
    if (nearBoundaryPoints_n[i] != n)
    {
        nearBoundaryPoints_n.insert(nearBoundaryPoints_n.begin() + i, n);
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::addNearBoundaryPoint_all(int n)
{

    if (nearBoundaryPoints_all.size() == 0)
    {
        nearBoundaryPoints_all.push_back(n);
        return;
    }
    int i = 0;
    for (i = 0; i < nearBoundaryPoints_all.size(); ++i)
    {
        if (n >= nearBoundaryPoints_all[i])
        {
            break;
        }
    }
    if (nearBoundaryPoints_all[i] != n)
    {
        nearBoundaryPoints_all.insert(nearBoundaryPoints_all.begin() + i, n);
    }
};

template <class PointType>
std::vector<int> PoissonSolver3d<PointType>::getNearBoundarypoints_n()
{
    return nearBoundaryPoints_n;
};

template <class PointType>
std::vector<int> PoissonSolver3d<PointType>::getNearBoundarypoints_all()
{
    return nearBoundaryPoints_all;
};

template <class PointType>
void PoissonSolver3d<PointType>::addBoundaryPoint(BoundaryPoint3d<PointType> bpoint,
                                                  DGeo::Point<PointType>     ipoint)
{

    int i = 0;

    if (boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].size() == 0)
    {
        boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].push_back(bpoint);
        return;
    }

    for (i = 0; i < boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].size(); ++i)
    {
        if (bpoint.pnum <= boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].pnum)
        {
            break;
        }
    }

    if (bpoint.pnum == boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].pnum)
    {
        if ((std::abs(bpoint.ipoint.x -
                 boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].ipoint.x) > eps_p) ||
            (std::abs(bpoint.ipoint.y -
                 boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].ipoint.y) > eps_p))
        {
            boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].edge.point1 = ipoint;
            boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].edge.point2 =
                boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].ipoint;
        }
        return;
    }

    if (bpoint.pnum != boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].pnum)
    {
        boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].insert(
            boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].begin() + i, bpoint);
        return;
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::addBoundaryPoint(BoundaryPoint3d<PointType> bpoint)
{

    int i = 0;

    if (boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].size() == 0)
    {
        boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].push_back(bpoint);
        return;
    }

    for (i = 0; i < boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].size(); ++i)
    {
        if (bpoint.pnum <= boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].pnum)
        {
            break;
        }
    }

    if (bpoint.pnum != boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1][i].pnum)
    {
        boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].insert(
            boundaryPoints[bpoint.coord.z - 1][bpoint.coord.y - 1].begin() + i, bpoint);
        return;
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::noteBoundaryOnMesh_3d(
    const std::shared_ptr<MeshContainer3d<PointType>>&            mesh,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
    const std::shared_ptr<BoundaryConditions>&                    boundaryConditions)
{

    int                        stepy, stepx;
    BoundaryPoint3d<PointType> bpoint, bpoint2;
    DGeo::Point<int>           mpoint1, mpoint2;
    DGeo::Point<PointType>     point1, point2;

    for (int z = 1; z < mesh->templNumb.GetNz() - 1; ++z)
    {

        int i = 1;
        int j = 1;

        while (mesh->templNumb(i, j, z) == -1)
        {
            if (j + 1 > mesh->templNumb.GetNcol() - 1)
            {
                ++i;
                j = 1;
            }
            else
            {
                ++j;
            }
        }

        searchBoundaryPointsIntersection_3d(mesh, boundaries, boundaryConditions,
                                            mesh->serialMeshData[mesh->templNumb(i, j, z)].z);

        // for z = z0
        for (int k = 0; k < boundaryPointsIntersectionCondition.size(); ++k)
        {
            point1    = boundaryPointsIntersectionCondition[k].point;
            mpoint1   = closestVertex_3d(mesh, point1);
            mpoint1.z = z; //костыль
            point2    = mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)];
            // steps
            stepx = 0;
            stepy = 0;
            if (std::abs(point2.x - point1.x) > eps_p)
            {
                stepx = sign(point1.x - point2.x);
            }
            if (std::abs(point2.y - point1.y) > eps_p)
            {
                stepy = sign(point1.y - point2.y);
            }
            if (stepx != 0)
            {
                mpoint2 = mesh->doStepX(mpoint1, stepx);
            }
            if (stepy != 0)
            {
                mpoint2 = mesh->doStepY(mpoint1, stepy);
            }

            if ((stepy != 0) && (stepx != 0))
            {
                double flags = 0;
            }

            if ((stepy == 0) && (stepx == 0))
            {
                mpoint2 = mpoint1;
            }
            int points = 0;
            if ((mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut == true) &&
                (mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)].isOut == true))
            {
                DGeo::Edge<PointType>  step_edge;
                DGeo::Point<PointType> point;
                step_edge.point1 = mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)];
                step_edge.point2 = mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)];
                points           = 0;
                for (int j = 0; j < boundaryPointsIntersectionDefault.size(); j++)
                {
                    point = boundaryPointsIntersectionDefault[j].point;

                    // is the point lay on the edge point1-point2
                    if (step_edge.IsOnEdge(point, eps_p) != 2)
                    {
                        points++;
                    }
                }
            }
            if ((points % 2 == 0) || (points == 1))
            {
                if (mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut == true)
                {
                    bpoint.edge  = boundaryPointsIntersectionCondition[k].edge;
                    bpoint.type  = boundaryPointsIntersectionCondition[k].type;
                    bpoint.value = boundaryPointsIntersectionCondition[k].value;
                    bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
                    bpoint.pnum  = mesh->templNumb(mpoint1.x, mpoint1.y, z);
                    bpoint.coord = mpoint1;
                }
                else if (mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)].isOut ==
                         true)
                {
                    bpoint.edge  = boundaryPointsIntersectionCondition[k].edge;
                    bpoint.type  = boundaryPointsIntersectionCondition[k].type;
                    bpoint.value = boundaryPointsIntersectionCondition[k].value;
                    bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
                    bpoint.pnum  = mesh->templNumb(mpoint2.x, mpoint2.y, z);
                    bpoint.coord = mpoint2;
                }
                else
                {
                    bpoint.edge  = boundaryPointsIntersectionCondition[k].edge;
                    bpoint.type  = boundaryPointsIntersectionCondition[k].type;
                    bpoint.value = boundaryPointsIntersectionCondition[k].value;
                    bpoint.pnum  = mesh->templNumb(mpoint1.x, mpoint1.y, z);
                    bpoint.coord = mpoint1;
                    bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
                    mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut = true;
                }
                addBoundaryPoint(bpoint, boundaryPointsIntersectionCondition[k].point);
            }
        }

        for (int k = 0; k < boundaryPointsIntersectionDefault.size(); ++k)
        {
            point1  = boundaryPointsIntersectionDefault[k].point;
            mpoint1 = closestVertex_3d(mesh, point1);
            point2  = mesh->mesh[z].meshData[mpoint1.y][mpoint1.x];
            // steps
            stepx = 0;
            stepy = 0;
            if (std::abs(point2.x - point1.x) > eps_p)
            {
                stepx = sign(point1.x - point2.x);
            }
            if (std::abs(point2.y - point1.y) > eps_p)
            {
                stepy = sign(point1.y - point2.y);
            }
            if (stepx != 0)
            {
                mpoint2 = mesh->doStepX(mpoint1, stepx);
            }
            if (stepy != 0)
            {
                mpoint2 = mesh->doStepY(mpoint1, stepy);
            }

            if ((stepy != 0) && (stepx != 0))
            {
                double flags = 0;
            }

            if ((stepy == 0) && (stepx == 0))
            {
                mpoint2 = mpoint1;
            }
            int points = 0;
            if ((mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut == true) &&
                (mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)].isOut == true))
            {
                DGeo::Edge<PointType>  step_edge;
                DGeo::Point<PointType> point;
                step_edge.point1 = mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)];
                step_edge.point2 = mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)];
                points           = 0;
                for (int j = 0; j < boundaryPointsIntersectionDefault.size(); j++)
                {
                    point = boundaryPointsIntersectionDefault[j].point;

                    // is the point lay on the edge point1-point2
                    if (step_edge.IsOnEdge(point, eps_p) != 2)
                    {
                        points++;
                    }
                }
            }
            if ((points % 2 == 0) || (points == 1))
            {
                if (mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut == true)
                {
                    bpoint.edge  = boundaryPointsIntersectionDefault[k].edge;
                    bpoint.type  = boundaryPointsIntersectionDefault[k].type;
                    bpoint.value = boundaryPointsIntersectionDefault[k].value;
                    bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
                    bpoint.pnum  = mesh->templNumb(mpoint1.x, mpoint1.y, z);
                    bpoint.coord = mpoint1;
                }
                else if (mesh->serialMeshData[mesh->templNumb(mpoint2.x, mpoint2.y, z)].isOut ==
                         true)
                {
                    bpoint.edge  = boundaryPointsIntersectionDefault[k].edge;
                    bpoint.type  = boundaryPointsIntersectionDefault[k].type;
                    bpoint.value = boundaryPointsIntersectionDefault[k].value;
                    bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
                    bpoint.pnum  = mesh->templNumb(mpoint2.x, mpoint2.y, z);
                    bpoint.coord = mpoint2;
                }
                else
                {
                    bpoint.edge  = boundaryPointsIntersectionDefault[k].edge;
                    bpoint.type  = boundaryPointsIntersectionDefault[k].type;
                    bpoint.value = boundaryPointsIntersectionDefault[k].value;
                    bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
                    bpoint.pnum  = mesh->mesh[z].meshData[mpoint1.y][mpoint1.x].Number;
                    bpoint.coord = mpoint1;
                    mesh->serialMeshData[mesh->templNumb(mpoint1.x, mpoint1.y, z)].isOut = true;
                }
                addBoundaryPoint(bpoint, boundaryPointsIntersectionDefault[k].point);
            }
        }

        for (int j = 1; j < mesh->templNumb.GetNrow() - 1; ++j)
        {
            for (int i = 1; i < mesh->templNumb.GetNcol() - 1; ++i)
            {
                if (mesh->templNumb(i, j, z) != -1)
                {
                    if (mesh->serialMeshData[mesh->templNumb(i, j, z)].isOut == true)
                    {

                        mpoint1.x      = i;
                        mpoint1.y      = j;
                        mpoint1.z      = z;
                        mpoint1.Number = mesh->templNumb(i, j, z);
                        bpoint2        = findBoundaryPoint(mpoint1);
                        if (bpoint2.type == 'f')
                        {
                            bpoint.type  = 'n';
                            bpoint.value = 0;
                            bpoint.bcnum = -1;
                            bpoint.pnum  = mesh->templNumb(i, j, z);
                            bpoint.coord = mpoint1;
                            addBoundaryPoint(bpoint);
                        }

                        /*
                        if (bpoint2.type != 'f'){
                        bpoint.edge = bpoint2.edge;
                        bpoint.type = bpoint2.type;
                        bpoint.value = bpoint2.value;
                        bpoint.bcnum = bpoint2.bcnum;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        else{
                        mpoint2 = mesh->doStepX(mpoint1, 1);
                        bpoint2 = findBoundaryPoint(mpoint2);
                        if (bpoint2.type != 'f'){
                        bpoint.edge = bpoint2.edge;
                        bpoint.type = bpoint2.type;
                        bpoint.value = bpoint2.value;
                        bpoint.bcnum = bpoint2.bcnum;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        else{
                        mpoint2 = mesh->doStepX(mpoint1, -1);
                        bpoint2 = findBoundaryPoint(mpoint2);
                        if (bpoint2.type != 'f'){
                        bpoint.edge = bpoint2.edge;
                        bpoint.type = bpoint2.type;
                        bpoint.value = bpoint2.value;
                        bpoint.bcnum = bpoint2.bcnum;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        else{
                        mpoint2 = mesh->doStepY(mpoint1, 1);
                        bpoint2 = findBoundaryPoint(mpoint2);
                        if (bpoint2.type != 'f'){
                        bpoint.edge = bpoint2.edge;
                        bpoint.type = bpoint2.type;
                        bpoint.value = bpoint2.value;
                        bpoint.bcnum = bpoint2.bcnum;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        else{
                        mpoint2 = mesh->doStepY(mpoint1, -1);
                        bpoint2 = findBoundaryPoint(mpoint2);
                        if (bpoint2.type != 'f'){
                        bpoint.edge = bpoint2.edge;
                        bpoint.type = bpoint2.type;
                        bpoint.value = bpoint2.value;
                        bpoint.bcnum = bpoint2.bcnum;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        else{
                        bpoint.type = 'n';
                        bpoint.value = 0;
                        bpoint.bcnum = -1;
                        bpoint.pnum = mesh->templNumb(i, j, z);
                        bpoint.coord = mpoint1;
                        }
                        }
                        }
                        }
                        }
                        addBoundaryPoint(bpoint);
                        */
                    }
                }
            }
        }
    }
    // i need to sort boundary points
    // std::sort(boundaryPoints.begin(), boundaryPoints.end(), compNum <PointType>);
};

template <class PointType>
void PoissonSolver3d<PointType>::createMatrix3d(
    const std::shared_ptr<MeshContainer3d<PointType>>& mesh, coeff3d left, coeff3d right,
    coeff3d up, coeff3d down, coeff3d outward, coeff3d deep, coeff3d middle)
{

    int size   = mesh->serialMeshData.size();
    int m      = 0;
    systemSize = size;
    int step;

    // int tempj;
    linVect = new double[size];
    x       = new double[size];
    for (int j = 0; j < size; ++j)
    {
        linVect[j] = 0; // right hand side
        x[j]       = 1; // start solution for iteration method
    }
    /*
    for (int j = 0; j < boundaryPoints.size(); ++j){
    x[boundaryPoints[j].pnum] = boundaryPoints[j].value;
    }
    */
    // loop for all points, building linear system
    DGeo::Point<int> tmp1, tmp2, tmp3, tmp4;
    int              row;
    double           h;
    int              flag = 0;

    double h_up, h_down, h_left, h_right, h_deep, h_outward;
    int    flag_h = 0;
    double h_bp_value, bp_bcnum;

    nonZeroLinVect.clear();
    NonZeroLinVect<PointType> tmp_nonZero;

    double                     x_, y_, hx_, hy_, z_, hz_;
    DGeo::Point<PointType>     norm;
    BoundaryPoint3d<PointType> bp;
    int                        c_i = -1; // current index

    for (int k = 0; k < mesh->templNumb.GetNz() - 1; k++)
    {
        tmp1.z = k;
        for (int i = 0; i < mesh->templNumb.GetNrow() - 1; i++)
        {
            tmp1.x = i;
            for (int j = 0; j < mesh->templNumb.GetNcol() - 1; j++)
            {
                tmp1.y = j;
                if (mesh->templNumb(i, j, k) != -1)
                {
                    tmp1.Number = mesh->templNumb(i, j, k);
                    if (mesh->serialMeshData[mesh->templNumb(i, j, k)].isOut == false)
                    {
                        row = mesh->templNumb(tmp1.x, tmp1.y, tmp1.z);
                        middles.push_back(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                        c_i++;

                        flag   = 0;
                        flag_h = 0;
                        // find distance to edges
                        // left
                        tmp2 = mesh->doStepX(tmp1, -1);
                        if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)].isOut ==
                            true)
                        {
                            bp     = findBoundaryPoint(tmp2);
                            h_left = distToEdgeX(
                                bp.edge,
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                            if (h_left < eps_h)
                            {
                                flag_h     = 1;
                                h_bp_value = bp.value;
                                bp_bcnum   = bp.bcnum;
                            }
                            if (bp.type == 'n')
                            {
                                h_left = distP2PX(
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                            }
                        }
                        else
                        {
                            h_left = distP2PX(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }
                        // right
                        tmp2 = mesh->doStepX(tmp1, 1);
                        if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)].isOut ==
                            true)
                        {
                            bp      = findBoundaryPoint(tmp2);
                            h_right = distToEdgeX(
                                bp.edge,
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                            if (h_right < eps_h)
                            {
                                flag_h     = 1;
                                h_bp_value = bp.value;
                                bp_bcnum   = bp.bcnum;
                            }
                            if (bp.type == 'n')
                            {
                                h_right = distP2PX(
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                            }
                        }
                        else
                        {
                            h_right = distP2PX(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }
                        // up
                        tmp2 = mesh->doStepY(tmp1, 1);
                        if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)].isOut ==
                            true)
                        {
                            bp   = findBoundaryPoint(tmp2);
                            h_up = distToEdgeY(
                                bp.edge,
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                            if (h_up < eps_h)
                            {
                                flag_h     = 1;
                                h_bp_value = bp.value;
                                bp_bcnum   = bp.bcnum;
                            }
                            if (bp.type == 'n')
                            {
                                h_up = distP2PY(
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                            }
                        }
                        else
                        {
                            h_up = distP2PY(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }

                        if (h_up != h_up)
                        {
                            h_up = distP2PY(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }
                        // down
                        tmp2 = mesh->doStepY(tmp1, -1);
                        if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)].isOut ==
                            true)
                        {
                            bp     = findBoundaryPoint(tmp2);
                            h_down = distToEdgeY(
                                bp.edge,
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                            if (h_down < eps_h)
                            {
                                flag_h     = 1;
                                h_bp_value = bp.value;
                                bp_bcnum   = bp.bcnum;
                            }
                            if (bp.type == 'n')
                            {
                                h_down = distP2PY(
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                            }
                        }
                        else
                        {
                            h_down = distP2PY(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }
                        if (h_down != h_down)
                        {
                            h_down = distP2PY(
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        }

                        // outward step = 1
                        tmp2 = mesh->doStepZ(tmp1, 1);
                        /*if (mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].isOut == true){
                        bp = findBoundaryPoint(mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].Number);
                        h_outward = distToEdgeZ(bp.edge,
                        mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x]);
                        if (h_outward < eps_h){
                        flag_h = 1;
                        h_bp_value = bp.value;
                        bp_bcnum = bp.bcnum;
                        }
                        }
                        else {*/
                        h_outward =
                            distP2PZ(mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                     mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        //}
                        // deep step = -1
                        tmp2 = mesh->doStepZ(tmp1, -1);
                        /*if (mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].isOut == true){
                        bp = findBoundaryPoint(mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].Number);
                        h_deep = distToEdgeZ(bp.edge, mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x]);
                        if (h_deep < eps_h){
                        flag_h = 1;
                        h_bp_value = bp.value;
                        bp_bcnum = bp.bcnum;
                        }
                        }
                        else {*/
                        h_deep =
                            distP2PZ(mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                     mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]);
                        //}
                        if ((h_up != h_up) || (h_down != h_down) || (h_right != h_right) ||
                            (h_left != h_left))
                        {
                            int stop = 1;
                        }

                        if (flag_h == 1)
                        {
                            // sys.setElement(row, row, 1);
                            c_middles.push_back(1);
                            linVect[row]      = linVect[row] + h_bp_value;
                            tmp_nonZero.bcnum = bp.bcnum;
                            tmp_nonZero.coef  = -1;
                            tmp_nonZero.row   = row;
                            nonZeroLinVect.push_back(tmp_nonZero);

                            addNearBoundaryPoint_all(mesh->templNumb(i, j, k));
                            c_rights.push_back(0);
                            c_lefts.push_back(0);
                            c_ups.push_back(0);
                            c_downs.push_back(0);
                            c_outwards.push_back(0);
                            c_deeps.push_back(0);

                            rights.push_back(0);
                            lefts.push_back(0);
                            ups.push_back(0);
                            downs.push_back(0);
                            deeps.push_back(0);
                            outwards.push_back(0);
                        }
                        else
                        {
                            x_ = mesh->serialMeshData[mesh->templNumb(i, j, k)].x;
                            y_ = mesh->serialMeshData[mesh->templNumb(i, j, k)].y;
                            z_ = mesh->serialMeshData[mesh->templNumb(i, j, k)].z;
                            // sys.setElement(row, row, middle(h_left, h_right, h_down, h_up, x_,
                            // y_));
                            c_middles.push_back(middle(h_left, h_right, h_down, h_up, h_deep,
                                                       h_outward, x_, y_, z_));
                            // right
                            step = 1;
                            h    = h_right;
                            tmp2 = mesh->doStepX(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            rights.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                if (bp.type == 'n')
                                {
                                    /*norm = normalVectToEdge(bp.edge);
                                    stepy = sign(norm.y);
                                    tmp3 = mesh->doStepY(tmp2, stepy);
                                    tmp4 = mesh->doStepY(tmp1, stepy);
                                    a = norm.x;
                                    b = norm.y;
                                    hx1 = distToEdgeX(bp.edge,mesh->meshData[tmp1.y][tmp1.x]);
                                    col2 = numbers[mesh->meshData[tmp4.y][tmp4.x].Number];
                                    //col3 = numbers[mesh->meshData[tmp3.y][tmp3.x].Number];
                                    */
                                    // 0 tmp1 row_row			+right(hx1,hy,x_,y_)
                                    // 1 tmp2 is out				linVect +bp.value()
                                    // 5 tmp3 is out/row_col3
                                    // 2 tmp4 row_col2
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // right(h_left, h_right,
                                    // h_down, h_up, x_, y_));
                                    c_middles[c_i] =
                                        c_middles[c_i] + right(h_left, h_right, h_down, h_up,
                                                               h_deep, h_outward, x_, y_, z_);
                                    c_rights.push_back(0);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       h * bp.value * step *
                                                           right(h_left, h_right, h_down, h_up,
                                                                 h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef =
                                        h * step * right(h_left, h_right, h_down, h_up, h_deep,
                                                         h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            bp.value * right(h_left, h_right, h_down, h_up, h_deep,
                                                             h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = right(h_left, h_right, h_down, h_up, h_deep,
                                                             h_outward, x_, y_, z_);
                                    tmp_nonZero.z   = tmp1.z;
                                    tmp_nonZero.row = row;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_rights.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, right(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_rights.push_back(right(h_left, h_right, h_down, h_up, h_deep,
                                                         h_outward, x_, y_, z_));
                            }
                            // left
                            step = -1;
                            h    = h_left;
                            tmp2 = mesh->doStepX(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            lefts.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                if (bp.type == 'n')
                                {
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // left(h_left, h_right, h_down,
                                    // h_up, x_, y_));
                                    c_lefts.push_back(0);
                                    c_middles[c_i] =
                                        c_middles[c_i] + left(h_left, h_right, h_down, h_up, h_deep,
                                                              h_outward, x_, y_, z_);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       h * bp.value * step *
                                                           left(h_left, h_right, h_down, h_up,
                                                                h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef =
                                        h * step * left(h_left, h_right, h_down, h_up, h_deep,
                                                        h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }

                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            bp.value * left(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = left(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_lefts.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, left(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_lefts.push_back(left(h_left, h_right, h_down, h_up, h_deep,
                                                       h_outward, x_, y_, z_));
                            }
                            // up
                            h    = h_up;
                            step = 1;
                            tmp2 = mesh->doStepY(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            ups.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                if (bp.type == 'n')
                                {
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // up(h_left, h_right, h_down,
                                    // h_up, x_, y_));
                                    c_middles[c_i] =
                                        c_middles[c_i] + up(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    c_ups.push_back(0);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            h * bp.value * step * up(h_left, h_right, h_down, h_up,
                                                                     h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = h * step * up(h_left, h_right, h_down, h_up,
                                                                     h_deep, h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       bp.value * up(h_left, h_right, h_down, h_up,
                                                                     h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = up(h_left, h_right, h_down, h_up, h_outward,
                                                          h_deep, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_ups.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, up(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_ups.push_back(up(h_left, h_right, h_down, h_up, h_deep, h_outward,
                                                   x_, y_, z_));
                            }
                            // down
                            h    = h_down;
                            step = -1;
                            tmp2 = mesh->doStepY(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            downs.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                // c_downs.push_back(0);
                                if (bp.type == 'n')
                                {
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // down(h_left, h_right, h_down,
                                    // h_up, x_, y_));
                                    c_middles[c_i] =
                                        c_middles[c_i] + down(h_left, h_right, h_down, h_up, h_deep,
                                                              h_outward, x_, y_, z_);
                                    c_downs.push_back(0);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       h * bp.value * step *
                                                           down(h_left, h_right, h_down, h_up,
                                                                h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef =
                                        h * step * down(h_left, h_right, h_down, h_up, h_deep,
                                                        h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            bp.value * down(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = down(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_downs.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, down(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_downs.push_back(down(h_left, h_right, h_down, h_up, h_deep,
                                                       h_outward, x_, y_, z_));
                            }
                            // outward
                            h    = h_outward;
                            step = 1;
                            tmp2 = mesh->doStepZ(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            outwards.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                // c_downs.push_back(0);
                                if (bp.type == 'n')
                                {
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // down(h_left, h_right, h_down,
                                    // h_up, x_, y_));
                                    c_middles[c_i] =
                                        c_middles[c_i] + outward(h_left, h_right, h_down, h_up,
                                                                 h_deep, h_outward, x_, y_, z_);
                                    c_outwards.push_back(0);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       h * bp.value * step *
                                                           outward(h_left, h_right, h_down, h_up,
                                                                   h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef =
                                        h * step * outward(h_left, h_right, h_down, h_up, h_deep,
                                                           h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            bp.value * outward(h_left, h_right, h_down, h_up,
                                                               h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = outward(h_left, h_right, h_down, h_up,
                                                               h_deep, h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_outwards.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, down(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_outwards.push_back(outward(h_left, h_right, h_down, h_up, h_deep,
                                                             h_outward, x_, y_, z_));
                            }
                            // deep
                            h    = h_down;
                            step = -1;
                            tmp2 = mesh->doStepZ(tmp1, step);
                            // col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
                            deeps.push_back(mesh->templNumb(tmp2.x, tmp2.y, tmp2.z));
                            if (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                    .isOut == true)
                            {
                                flag = 1;
                                bp   = findBoundaryPoint(tmp2);
                                // c_downs.push_back(0);
                                if (bp.type == 'n')
                                {
                                    // sys.setElement(row, row, sys.getElement(row, row) +
                                    // down(h_left, h_right, h_down,
                                    // h_up, x_, y_));
                                    c_middles[c_i] =
                                        c_middles[c_i] + deep(h_left, h_right, h_down, h_up, h_deep,
                                                              h_outward, x_, y_, z_);
                                    c_deeps.push_back(0);
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] = linVect[row] -
                                                       h * bp.value * step *
                                                           deep(h_left, h_right, h_down, h_up,
                                                                h_deep, h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef =
                                        h * step * deep(h_left, h_right, h_down, h_up, h_deep,
                                                        h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    addNearBoundaryPoint_n(mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                                if (bp.type == 'd')
                                {
                                    if (std::abs(bp.value) > 0.00001)
                                    {
                                        linVect[row] =
                                            linVect[row] -
                                            bp.value * deep(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    }
                                    tmp_nonZero.bcnum = bp.bcnum;
                                    tmp_nonZero.coef  = deep(h_left, h_right, h_down, h_up, h_deep,
                                                            h_outward, x_, y_, z_);
                                    tmp_nonZero.row = row;
                                    tmp_nonZero.z   = tmp1.z;
                                    nonZeroLinVect.push_back(tmp_nonZero);

                                    c_deeps.push_back(0);
                                    addNearBoundaryPoint_all(
                                        mesh->templNumb(tmp1.x, tmp1.y, tmp1.z));
                                }
                            }
                            else
                            {
                                // sys.setElement(row, col, down(h_left, h_right, h_down, h_up, x_,
                                // y_));
                                c_deeps.push_back(deep(h_left, h_right, h_down, h_up, h_deep,
                                                       h_outward, x_, y_, z_));
                            }
                        }
                    }
                }
            }
        }
    }
    std::sort(nearBoundaryPoints_all.begin(), nearBoundaryPoints_all.end());
    std::sort(nearBoundaryPoints_n.begin(), nearBoundaryPoints_n.end());
};

template <class PointType>
void PoissonSolver3d<PointType>::InitSolver(
    const std::shared_ptr<GridData3d<PointType>>&                 gridData,
    const std::shared_ptr<MeshContainer3d<PointType>>&            mesh,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
    const std::shared_ptr<BoundaryConditions>&                    boundaryConditions)
{

    boundaryPoints.clear();
    nearBoundaryPoints_all.clear();
    nearBoundaryPoints_n.clear();

    middles.clear();
    ups.clear();
    downs.clear();
    c_rights.clear();
    c_middles.clear();
    lefts.clear();
    rights.clear();
    c_downs.clear();
    c_lefts.clear();
    c_ups.clear();
    c_outwards.clear();
    c_deeps.clear();
    outwards.clear();
    deeps.clear();
    nSliceBPoints.clear();
    boundaryPointsIntersectionCondition.clear();
    boundaryPointsIntersectionDefault.clear();

    nSlicePoints = mesh->serialMeshData.size() / mesh->mesh.size();

    noteBoundaryOnMesh_3d(mesh, boundaries, boundaryConditions);

    int k = 0;

    for (int i = 0; i < mesh->serialMeshData.size(); ++i)
    {
        mesh->serialMeshData[i].Number = i;
    }
    polar = false;

    BoundaryPoint3d<PointType> tmp_bp;
    tmp_bp.type  = 'n';
    tmp_bp.value = 0;
    tmp_bp.bcnum = -1;
    // note front and back sides of mesh like neumann boundary
    for (int i = 0; i < mesh->templNumb.GetNrow() - 1; i++)
    {
        for (int j = 0; j < mesh->templNumb.GetNcol() - 1; j++)
        {
            if (mesh->templNumb(i, j, 1) != -1)
            {
                mesh->serialMeshData[mesh->templNumb(i, j, 1)].isOut = true;
            }
            if (mesh->templNumb(i, j, mesh->templNumb.GetNz() - 2) != -1)
            {
                mesh->serialMeshData[mesh->templNumb(i, j, mesh->templNumb.GetNz() - 2)].isOut =
                    true;
            }
        }
    }

    std::vector<std::vector<BoundaryPoint3d<PointType>>> ttmp_v;
    ttmp_v.resize(mesh->templNumb.GetNcol() - 2);
    boundaryPoints.resize(mesh->templNumb.GetNz() - 2, ttmp_v);

    double x, y;
    double r = 0.015;

    for (int k = 1; k < mesh->templNumb.GetNz() - 1; ++k)
    {
        int i = 1;
        int j = 1;
        while (mesh->templNumb(i, j, k) == -1)
        {
            if (j + 1 > mesh->templNumb.GetNcol() - 1)
            {
                ++i;
                j = 1;
            }
            else
            {
                ++j;
            }
        }

        if (std::abs(mesh->serialMeshData[mesh->templNumb(i, j, k)].z - 0.0) < eps_h)
        {
            for (int i = 1; i < mesh->templNumb.GetNrow() - 1; i++)
            {
                for (int j = 1; j < mesh->templNumb.GetNcol() - 1; j++)
                {
                    if (mesh->templNumb(i, j, k) != -1)
                    {
                        x = mesh->serialMeshData[mesh->templNumb(i, j, k)].x;
                        y = mesh->serialMeshData[mesh->templNumb(i, j, k)].y;
                        if (sqrt(x * x + y * y) > r)
                        {
                            tmp_bp.bcnum   = 3;
                            tmp_bp.pnum    = mesh->templNumb(i, j, k);
                            tmp_bp.type    = 'd';
                            tmp_bp.value   = 30;
                            tmp_bp.coord.x = i;
                            tmp_bp.coord.y = j;
                            tmp_bp.coord.z = k;
                            addBoundaryPoint(tmp_bp);
                            mesh->serialMeshData[mesh->templNumb(i, j, k)].isOut = true;
                        }
                    }
                }
            }
        }
    }

    createMatrix3d(mesh, left3d, right3d, up3d, down3d, outward3d, deep3d, middle3d);

    std::vector<PointType> u;
    u.resize(systemSize, 0);
    U.push_back(u); // charge
    solve_offset(u, mesh, boundaryConditions);
    U.push_back(u); // offset
    std::vector<std::vector<double>> zArray;
    std::vector<PropertyCondition>   local_list;
    local_list = boundaryConditions->GetPropertyConditionsList();
    double amplitude;
    int    status;

    for (int i = 0; i < boundaryConditions->PropertyConditionListSize(); ++i)
    {
        zArray = local_list[i].Get_zArray();
        for (int j = 0; j < zArray.size(); ++j)
        {
            amplitude = boundaryConditions->GetPotentialAmplitude(
                i, 0.5 * zArray[j][0] + 0.5 * zArray[j][1], status);
            if (std::abs(amplitude) > eps_tolerance)
            {
                u.resize(systemSize);
                solve_one_rotating(gridData->Getrho(), u, mesh, boundaryConditions, i, zArray[j][0],
                                   zArray[j][1]);
                U.push_back(u);
                u.clear();
            }
        }
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::FieldSimulate(
    const std::shared_ptr<GridData3d<PointType>>&      gridData,
    const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
    const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t)
{

    // solve_charge(gridData->Getrho(), U[0], mesh, boundaryConditions);

    sum_solutions(U, gridData->GetV(), boundaryConditions, t);
    extrapolate_potential(gridData->GetV(), mesh);

    fieldCalculation_Simple_FDE_3d(gridData->GetV(), mesh, 'x', gridData->Get_Ex());
    fieldCalculation_Simple_FDE_3d(gridData->GetV(), mesh, 'y', gridData->Get_Ey());
    fieldCalculation_Simple_FDE_3d(gridData->GetV(), mesh, 'z', gridData->Get_Ez());
};

template <class PointType>
void PoissonSolver3d<PointType>::sum_solutions(
    std::vector<std::vector<PointType>>& U, std::vector<PointType>& V,
    const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t)
{

    std::vector<PointType>           koefs(U.size() - 2);
    std::vector<std::vector<double>> zArray;
    int                              status;
    PointType                        koef;
    std::vector<PropertyCondition>   local_list;
    local_list = boundaryConditions->GetPropertyConditionsList();
    double amplitude;

    for (int i = 0; i < boundaryConditions->PropertyConditionListSize(); ++i)
    {
        zArray = local_list[i].Get_zArray();
        for (int j = 0; j < zArray.size(); ++j)
        {
            amplitude = boundaryConditions->GetPotentialAmplitude(
                i, 0.5 * zArray[j][0] + 0.5 * zArray[j][1], status);
            if (std::abs(amplitude) > eps_tolerance)
            {
                koef = boundaryConditions->GetPotential(
                           i, t, 0.5 * zArray[j][0] + 0.5 * zArray[j][1], status) /
                       boundaryConditions->GetPotentialAmplitude(
                           i, 0.5 * zArray[j][0] + 0.5 * zArray[j][1], status);
                koefs.push_back(koef);
            }
        }
    }

    for (int j = 0; j < V.size(); ++j)
    {
        V[j] = U[0][j] + U[1][j];
    }

    for (int i = 2; i < U.size(); ++i)
    {
        for (int j = 0; j < V.size(); ++j)
        {
            V[j] = V[j] + U[i][j] * koefs[i - 2];
        }
    }
}

template <class PointType>
void PoissonSolver3d<PointType>::solve_one_rotating(
    std::vector<PointType>& rho, std::vector<PointType>& V,
    const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
    const std::shared_ptr<BoundaryConditions>& boundaryConditions, int cond_num, double z_begin,
    double z_end)
{

    double           pi   = 3.14159265358979323;
    int              size = systemSize;
    std::vector<int> strElements;
    double           tempValue;
    double           eps0     = commtools::VACUUM_PERMITTIVITY();
    int              flag     = 0;
    double           tmp_sum  = 0;
    double           diff_sum = 0;
    double*          vect     = new double[size];

    int status = 0;

    // get z coordinates
    std::vector<double> z_s;
    for (int k = 1; k < mesh->templNumb.GetNz() - 1; ++k)
    {
        int i = 1;
        int j = 1;
        while (mesh->templNumb(i, j, k) == -1)
        {
            if (j + 1 > mesh->templNumb.GetNcol() - 1)
            {
                ++i;
                j = 1;
            }
            else
            {
                ++j;
            }
        }
        z_s.push_back(mesh->serialMeshData[mesh->templNumb(i, j, k)].z);
    }

    std::vector<std::vector<PointType>> boundConds;
    std::vector<PointType>              tmp_boundCond;
    boundConds.clear();
    for (int j = 0; j < z_s.size(); ++j)
    {
        for (int i = 0; i < boundCondNum; ++i)
        {
            if ((i == cond_num) && (z_s[j] > z_begin - eps_p) && (z_s[j] < z_end + eps_p))
            {
                tmp_boundCond.push_back(
                    boundaryConditions->GetPotentialAmplitude(i, z_s[j], status));
            }
            else
            {
                tmp_boundCond.push_back(0);
            }
        }
        boundConds.push_back(tmp_boundCond);
        tmp_boundCond.clear();
    }

    for (int i = 0; i < size; ++i)
    {
        linVect[i] = 0;
    }

    for (int i = 0; i < nonZeroLinVect.size(); ++i)
    {
        if (nonZeroLinVect[i].bcnum != -1)
        {
            linVect[nonZeroLinVect[i].row] =
                linVect[nonZeroLinVect[i].row] -
                nonZeroLinVect[i].coef *
                    boundConds[nonZeroLinVect[i].z - 2][nonZeroLinVect[i].bcnum];
        }
    }

    double step_z = z_s[1] - z_s[0];

    for (int k = 0; k < boundaryPoints.size(); ++k)
    {
        for (int j = 0; j < boundaryPoints[k].size(); ++j)
        {
            for (int i = 0; i < boundaryPoints[k][j].size(); ++i)
            {
                x[boundaryPoints[k][j][i].pnum] =
                    boundConds[round(mesh->serialMeshData[boundaryPoints[k][j][i].pnum].z / step_z)]
                              [boundaryPoints[k][j][i].bcnum];
            }
        }
    }

    for (int i = 0; i < rho.size(); ++i)
    {
        vect[i] = linVect[i]; // -rho[i] / eps0;
    }

    std::vector<double> r, r_, v_, p, s, t_;
    for (int i = 0; i < size; ++i)
    {
        r.push_back(0);
        r_.push_back(0);
        v_.push_back(0);
        p.push_back(0);
        s.push_back(0);
        t_.push_back(0);
    }

    double ro, alpha, omega, beta, tmp_sum2;

    for (int i = 0; i < middles.size(); ++i)
    {
        v_[middles[i]] = 0;
        p[middles[i]]  = 0;
        r[middles[i]] =
            vect[middles[i]] -
            (x[middles[i]] * c_middles[i] + x[rights[i]] * c_rights[i] + x[lefts[i]] * c_lefts[i] +
             x[ups[i]] * c_ups[i] + x[downs[i]] * c_downs[i] + x[deeps[i]] * c_deeps[i] +
             x[outwards[i]] * c_outwards[i]);
        r_[middles[i]] = r[middles[i]];
    }
    ro    = 1;
    alpha = 1;
    omega = 1;
    while (flag == 0)
    {

        beta = alpha / (ro * omega);
        ro   = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            ro = ro + r_[middles[i]] * r[middles[i]];
        }

        beta    = beta * ro;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            p[middles[i]] = r[middles[i]] + beta * (p[middles[i]] - omega * v_[middles[i]]);
        }
        for (int i = 0; i < middles.size(); ++i)
        {
            v_[middles[i]] =
                (p[middles[i]] * c_middles[i] + p[rights[i]] * c_rights[i] +
                 p[lefts[i]] * c_lefts[i] + p[ups[i]] * c_ups[i] + p[downs[i]] * c_downs[i] +
                 p[deeps[i]] * c_deeps[i] + p[outwards[i]] * c_outwards[i]);
            tmp_sum = tmp_sum + r_[middles[i]] * v_[middles[i]];
        }

        alpha   = ro / tmp_sum;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            s[middles[i]] = r[middles[i]] - alpha * v_[middles[i]];
            tmp_sum       = tmp_sum + s[middles[i]];
        }

        if (std::abs(tmp_sum) < eps_tolerance)
        {
            flag = 1;
        }

        for (int i = 0; i < middles.size(); ++i)
        {
            t_[middles[i]] =
                (s[middles[i]] * c_middles[i] + s[rights[i]] * c_rights[i] +
                 s[lefts[i]] * c_lefts[i] + s[ups[i]] * c_ups[i] + s[downs[i]] * c_downs[i] +
                 s[deeps[i]] * c_deeps[i] + s[outwards[i]] * c_outwards[i]);
        }

        tmp_sum  = 0;
        tmp_sum2 = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            tmp_sum  = tmp_sum + t_[middles[i]] * s[middles[i]];
            tmp_sum2 = tmp_sum2 + t_[middles[i]] * t_[middles[i]];
        }
        omega = tmp_sum / tmp_sum2;

        for (int i = 0; i < middles.size(); ++i)
        {
            r[middles[i]] = s[middles[i]] - omega * t_[middles[i]];
        }
        diff_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            diff_sum      = diff_sum + std::abs(omega * s[middles[i]] + alpha * p[middles[i]]);
            x[middles[i]] = x[middles[i]] + omega * s[middles[i]] + alpha * p[middles[i]];
        }
        if (diff_sum < eps_tolerance)
        {
            flag = 1;
        }
    }

    /*
    int n = 0;
    double sum;
    while (flag == 0){
    diff_sum = 0;
    tmp_sum = 0;
    for (int k = 0; k < middles.size(); ++k){
    tmp_sum = 0;
    n = middles[k];
    tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] + x[ups[k]] * c_ups[k] +
    x[downs[k]] * c_downs[k]
    + x[deeps[k]] * c_deeps[k] + x[outwards[k]] * c_outwards[k])*w / c_middles[k];
    tmp_sum = tmp_sum + (1 - w)*x[n] + w*vect[n] / c_middles[k];
    if ((x[n]> eps_tolerance*0.1) || (x[n] < -eps_tolerance*0.1)){
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum) / std::abs(x[n]);
    }
    else {
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum);
    }
    x[n] = tmp_sum;
    }
    if (std::abs(diff_sum) < eps_tolerance) flag = 1;
    }
    */

    // write solution in gridData
    int m = 0;
    for (int i = 0; i < size; ++i)
    {
        V[i] = x[i];
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::solve_charge(
    std::vector<PointType>& rho, std::vector<PointType>& V,
    const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
    const std::shared_ptr<BoundaryConditions>&         boundaryConditions)
{

    double           pi   = 3.14159265358979323;
    int              size = systemSize;
    std::vector<int> strElements;
    double           tempValue;
    double           eps0     = commtools::VACUUM_PERMITTIVITY();
    int              flag     = 0;
    double           tmp_sum  = 0;
    double           diff_sum = 0;
    double*          vect     = new double[size];

    for (int i = 0; i < size; ++i)
    {
        linVect[i] = 0;
    }

    for (int k = 0; k < boundaryPoints.size(); ++k)
    {
        for (int j = 0; j < boundaryPoints[k].size(); ++j)
        {
            for (int i = 0; i < boundaryPoints[k][j].size(); ++i)
            {
                x[boundaryPoints[k][j][i].pnum] = 0;
            }
        }
    }

    for (int i = 0; i < rho.size(); ++i)
    {
        vect[i] = linVect[i] - rho[i] / eps0;
    }

    std::vector<double> r, r_, v_, p, s, t_;
    for (int i = 0; i < size; ++i)
    {
        r.push_back(0);
        r_.push_back(0);
        v_.push_back(0);
        p.push_back(0);
        s.push_back(0);
        t_.push_back(0);
    }

    double ro, alpha, omega, beta, tmp_sum2;

    for (int i = 0; i < middles.size(); ++i)
    {
        v_[middles[i]] = 0;
        p[middles[i]]  = 0;
        r[middles[i]] =
            vect[middles[i]] -
            (x[middles[i]] * c_middles[i] + x[rights[i]] * c_rights[i] + x[lefts[i]] * c_lefts[i] +
             x[ups[i]] * c_ups[i] + x[downs[i]] * c_downs[i] + x[deeps[i]] * c_deeps[i] +
             x[outwards[i]] * c_outwards[i]);
        r_[middles[i]] = r[middles[i]];
    }
    ro    = 1;
    alpha = 1;
    omega = 1;
    while (flag == 0)
    {

        beta = alpha / (ro * omega);
        ro   = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            ro = ro + r_[middles[i]] * r[middles[i]];
        }

        beta    = beta * ro;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            p[middles[i]] = r[middles[i]] + beta * (p[middles[i]] - omega * v_[middles[i]]);
        }
        for (int i = 0; i < middles.size(); ++i)
        {
            v_[middles[i]] =
                (p[middles[i]] * c_middles[i] + p[rights[i]] * c_rights[i] +
                 p[lefts[i]] * c_lefts[i] + p[ups[i]] * c_ups[i] + p[downs[i]] * c_downs[i] +
                 p[deeps[i]] * c_deeps[i] + p[outwards[i]] * c_outwards[i]);
            tmp_sum = tmp_sum + r_[middles[i]] * v_[middles[i]];
        }

        alpha   = ro / tmp_sum;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            s[middles[i]] = r[middles[i]] - alpha * v_[middles[i]];
            tmp_sum       = tmp_sum + s[middles[i]];
        }

        if (std::abs(tmp_sum) < eps_tolerance)
        {
            flag = 1;
        }

        for (int i = 0; i < middles.size(); ++i)
        {
            t_[middles[i]] =
                (s[middles[i]] * c_middles[i] + s[rights[i]] * c_rights[i] +
                 s[lefts[i]] * c_lefts[i] + s[ups[i]] * c_ups[i] + s[downs[i]] * c_downs[i] +
                 s[deeps[i]] * c_deeps[i] + s[outwards[i]] * c_outwards[i]);
        }

        tmp_sum  = 0;
        tmp_sum2 = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            tmp_sum  = tmp_sum + t_[middles[i]] * s[middles[i]];
            tmp_sum2 = tmp_sum2 + t_[middles[i]] * t_[middles[i]];
        }
        omega = tmp_sum / tmp_sum2;

        for (int i = 0; i < middles.size(); ++i)
        {
            r[middles[i]] = s[middles[i]] - omega * t_[middles[i]];
        }
        diff_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            diff_sum      = diff_sum + std::abs(omega * s[middles[i]] + alpha * p[middles[i]]);
            x[middles[i]] = x[middles[i]] + omega * s[middles[i]] + alpha * p[middles[i]];
        }
        if (diff_sum < eps_tolerance)
        {
            flag = 1;
        }
    }

    /*
    int n = 0;
    double sum;
    while (flag == 0){
    diff_sum = 0;
    tmp_sum = 0;
    for (int k = 0; k < middles.size(); ++k){
    tmp_sum = 0;
    n = middles[k];
    tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] + x[ups[k]] * c_ups[k] +
    x[downs[k]] * c_downs[k]
    + x[deeps[k]] * c_deeps[k] + x[outwards[k]] * c_outwards[k])*w / c_middles[k];
    tmp_sum = tmp_sum + (1 - w)*x[n] + w*vect[n] / c_middles[k];
    if ((x[n]> eps_tolerance*0.1) || (x[n] < -eps_tolerance*0.1)){
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum) / std::abs(x[n]);
    }
    else {
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum);
    }
    x[n] = tmp_sum;
    }
    if (std::abs(diff_sum) < eps_tolerance) flag = 1;
    }
    */

    // write solution in gridData
    int m = 0;
    for (int i = 0; i < size; ++i)
    {
        V[i] = x[i];
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::solve_offset(
    std::vector<PointType>& V, const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
    const std::shared_ptr<BoundaryConditions>& boundaryConditions)
{

    double           pi   = 3.14159265358979323;
    int              size = systemSize;
    std::vector<int> strElements;
    double           tempValue;
    double           eps0     = commtools::VACUUM_PERMITTIVITY();
    int              flag     = 0;
    double           tmp_sum  = 0;
    double           diff_sum = 0;
    double*          vect     = new double[size];

    int status = 0;

    // get z coordinates
    std::vector<double> z_s;
    for (int k = 1; k < mesh->templNumb.GetNz() - 1; ++k)
    {
        int i = 1;
        int j = 1;
        while (mesh->templNumb(i, j, k) == -1)
        {
            if (j + 1 > mesh->templNumb.GetNcol() - 1)
            {
                ++i;
                j = 1;
            }
            else
            {
                ++j;
            }
        }
        z_s.push_back(mesh->serialMeshData[mesh->templNumb(i, j, k)].z);
    }

    std::vector<std::vector<PointType>> boundConds;
    std::vector<PointType>              tmp_boundCond;
    boundConds.clear();
    for (int j = 0; j < z_s.size(); ++j)
    {
        for (int i = 0; i < boundCondNum; ++i)
        {
            tmp_boundCond.push_back(boundaryConditions->GetPotentialOffset(i, z_s[j], status));
        }
        boundConds.push_back(tmp_boundCond);
        tmp_boundCond.clear();
    }

    for (int i = 0; i < size; ++i)
    {
        linVect[i] = 0;
    }

    for (int i = 0; i < nonZeroLinVect.size(); ++i)
    {
        if (nonZeroLinVect[i].bcnum != -1)
        {
            linVect[nonZeroLinVect[i].row] =
                linVect[nonZeroLinVect[i].row] -
                nonZeroLinVect[i].coef *
                    boundConds[nonZeroLinVect[i].z - 2][nonZeroLinVect[i].bcnum];
        }
    }

    double step_z = z_s[1] - z_s[0];

    for (int k = 0; k < boundaryPoints.size(); ++k)
    {
        for (int j = 0; j < boundaryPoints[k].size(); ++j)
        {
            for (int i = 0; i < boundaryPoints[k][j].size(); ++i)
            {
                x[boundaryPoints[k][j][i].pnum] =
                    boundConds[round(mesh->serialMeshData[boundaryPoints[k][j][i].pnum].z / step_z)]
                              [boundaryPoints[k][j][i].bcnum];
            }
        }
    }

    for (int i = 0; i < V.size(); ++i)
    {
        vect[i] = linVect[i];
    }

    std::vector<double> r, r_, v_, p, s, t_;
    for (int i = 0; i < size; ++i)
    {
        r.push_back(0);
        r_.push_back(0);
        v_.push_back(0);
        p.push_back(0);
        s.push_back(0);
        t_.push_back(0);
    }

    double ro, alpha, omega, beta, tmp_sum2;

    for (int i = 0; i < middles.size(); ++i)
    {
        v_[middles[i]] = 0;
        p[middles[i]]  = 0;
        r[middles[i]] =
            vect[middles[i]] -
            (x[middles[i]] * c_middles[i] + x[rights[i]] * c_rights[i] + x[lefts[i]] * c_lefts[i] +
             x[ups[i]] * c_ups[i] + x[downs[i]] * c_downs[i] + x[deeps[i]] * c_deeps[i] +
             x[outwards[i]] * c_outwards[i]);
        r_[middles[i]] = r[middles[i]];
    }
    ro    = 1;
    alpha = 1;
    omega = 1;
    while (flag == 0)
    {

        beta = alpha / (ro * omega);
        ro   = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            ro = ro + r_[middles[i]] * r[middles[i]];
        }

        beta    = beta * ro;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            p[middles[i]] = r[middles[i]] + beta * (p[middles[i]] - omega * v_[middles[i]]);
        }
        for (int i = 0; i < middles.size(); ++i)
        {
            v_[middles[i]] =
                (p[middles[i]] * c_middles[i] + p[rights[i]] * c_rights[i] +
                 p[lefts[i]] * c_lefts[i] + p[ups[i]] * c_ups[i] + p[downs[i]] * c_downs[i] +
                 p[deeps[i]] * c_deeps[i] + p[outwards[i]] * c_outwards[i]);
            tmp_sum = tmp_sum + r_[middles[i]] * v_[middles[i]];
        }

        alpha   = ro / tmp_sum;
        tmp_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            s[middles[i]] = r[middles[i]] - alpha * v_[middles[i]];
            tmp_sum       = tmp_sum + s[middles[i]];
        }

        if (std::abs(tmp_sum) < eps_tolerance)
        {
            flag = 1;
        }

        for (int i = 0; i < middles.size(); ++i)
        {
            t_[middles[i]] =
                (s[middles[i]] * c_middles[i] + s[rights[i]] * c_rights[i] +
                 s[lefts[i]] * c_lefts[i] + s[ups[i]] * c_ups[i] + s[downs[i]] * c_downs[i] +
                 s[deeps[i]] * c_deeps[i] + s[outwards[i]] * c_outwards[i]);
        }

        tmp_sum  = 0;
        tmp_sum2 = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            tmp_sum  = tmp_sum + t_[middles[i]] * s[middles[i]];
            tmp_sum2 = tmp_sum2 + t_[middles[i]] * t_[middles[i]];
        }
        omega = tmp_sum / tmp_sum2;

        for (int i = 0; i < middles.size(); ++i)
        {
            r[middles[i]] = s[middles[i]] - omega * t_[middles[i]];
        }
        diff_sum = 0;
        for (int i = 0; i < middles.size(); ++i)
        {
            diff_sum      = diff_sum + std::abs(omega * s[middles[i]] + alpha * p[middles[i]]);
            x[middles[i]] = x[middles[i]] + omega * s[middles[i]] + alpha * p[middles[i]];
        }
        if (diff_sum < eps_tolerance)
        {
            flag = 1;
        }
    }

    /*
    int n = 0;
    double sum;
    while (flag == 0){
    diff_sum = 0;
    tmp_sum = 0;
    for (int k = 0; k < middles.size(); ++k){
    tmp_sum = 0;
    n = middles[k];
    tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] + x[ups[k]] * c_ups[k] +
    x[downs[k]] * c_downs[k]
    + x[deeps[k]] * c_deeps[k] + x[outwards[k]] * c_outwards[k])*w / c_middles[k];
    tmp_sum = tmp_sum + (1 - w)*x[n] + w*vect[n] / c_middles[k];
    if ((x[n]> eps_tolerance*0.1) || (x[n] < -eps_tolerance*0.1)){
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum) / std::abs(x[n]);
    }
    else {
    diff_sum = diff_sum + std::abs(x[n] - tmp_sum);
    }
    x[n] = tmp_sum;
    }
    if (std::abs(diff_sum) < eps_tolerance) flag = 1;
    }
    */

    // write solution in gridData
    int m = 0;
    for (int i = 0; i < size; ++i)
    {
        V[i] = x[i];
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::extrapolate_potential(
    std::vector<PointType>& V, const std::shared_ptr<MeshContainer3d<PointType>>& mesh)
{
    // extrapolation potential
    DGeo::Point<int>              tmp, tmp1, tmp2;
    PointType                     tmp_V, v, h_, v1, v2, h1, h2;
    BoundaryPoint3d<PointType>    bp;
    std::vector<DGeo::Point<int>> strange_points;
    double                        another_h;
    int                           nm;

    for (int k = 1; k < mesh->templNumb.GetNz() - 1; k++)
    {
        for (int i = 1; i < mesh->templNumb.GetNrow() - 1; i++)
        {
            for (int j = 1; j < mesh->templNumb.GetNcol() - 1; j++)
            {
                if (mesh->templNumb(i, j, k) != -1)
                {
                    if (mesh->serialMeshData[mesh->templNumb(i, j, k)].isOut == true)
                    {
                        tmp.x      = i;
                        tmp.y      = j;
                        tmp.z      = k;
                        tmp.Number = mesh->templNumb(i, j, k);
                        bp         = findBoundaryPoint(tmp);
                        nm         = 0;
                        tmp_V      = 0;

                        if (bp.type == 'n')
                        {
                            tmp1 = mesh->doStepX(tmp, 1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepX(tmp1, 1);
                                h1   = distP2PX(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                { //потому что нет резаных ячеек для неймана
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }
                            tmp1 = mesh->doStepX(tmp, -1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepX(tmp1, -1);
                                h1   = distP2PX(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                {
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }
                            tmp1 = mesh->doStepY(tmp, 1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepY(tmp1, 1);
                                h1   = distP2PY(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                {
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PY(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }
                            tmp1 = mesh->doStepY(tmp, -1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepY(tmp1, -1);
                                h1   = distP2PY(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                {
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PY(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }

                            tmp1 = mesh->doStepZ(tmp, 1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepZ(tmp1, 1);
                                h1   = distP2PZ(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                {
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PZ(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }
                            tmp1 = mesh->doStepZ(tmp, -1);
                            if (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                    .isOut == false)
                            {
                                tmp2 = mesh->doStepZ(tmp1, -1);
                                h1   = distP2PZ(
                                    mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                    mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]);
                                v1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z))
                                {
                                    v2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                    h2 = distP2PZ(mesh->serialMeshData[mesh->templNumb(
                                                      tmp1.x, tmp1.y, tmp1.z)],
                                                  mesh->serialMeshData[mesh->templNumb(
                                                      tmp2.x, tmp2.y, tmp2.z)]);
                                    tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
                                    nm++;
                                }
                            }

                            if (nm != 0)
                            {
                                V[mesh->templNumb(tmp.x, tmp.y, tmp.z)] = tmp_V / (1.0 * nm);
                            }
                            else
                            {
                                strange_points.push_back(tmp);
                            }
                        }
                        else
                        {
                            DGeo::Edge<PointType> points_edge;
                            DGeo::Edge<PointType> edge;
                            points_edge.point1 =
                                mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                            tmp1 = mesh->doStepX(tmp, 1);
                            points_edge.point2 =
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                            if ((mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                     .isOut == false))
                            {
                                if (isIntersection(bp.edge, points_edge, eps_p))
                                {
                                    h_ = distToEdgeX(bp.edge, mesh->serialMeshData[mesh->templNumb(
                                                                  tmp1.x, tmp1.y, tmp1.z)]);
                                    h1 = distP2PX(
                                        mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                        mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y,
                                                                             tmp1.z)]);
                                    v1 = V[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                                    v2 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                    v  = h1 * (v1 - v2) / h_ + v2;
                                    nm++;
                                    if (h_ < eps_h)
                                    {
                                        tmp2 = mesh->doStepX(tmp1, 1);
                                        v2   = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                        h2   = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                          tmp2.x, tmp2.y, tmp2.z)],
                                                      mesh->serialMeshData[mesh->templNumb(
                                                          tmp1.x, tmp1.y, tmp1.z)]);
                                        v = (h1 - h_) * (v1 - v2) / h2 + v1;
                                    }
                                    tmp_V = tmp_V + v;
                                }
                            }
                            tmp1 = mesh->doStepX(tmp, -1);
                            points_edge.point2 =
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                            if ((mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                     .isOut == false))
                            {
                                if (isIntersection(bp.edge, points_edge, eps_p))
                                {
                                    h_ = distToEdgeX(bp.edge, mesh->serialMeshData[mesh->templNumb(
                                                                  tmp1.x, tmp1.y, tmp1.z)]);
                                    h1 = distP2PX(
                                        mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                        mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y,
                                                                             tmp1.z)]);
                                    v1 = V[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                                    v2 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                    v  = h1 * (v1 - v2) / h_ + v2;
                                    nm++;
                                    if (h_ < eps_h)
                                    {
                                        tmp2 = mesh->doStepX(tmp1, -1);
                                        v2   = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                        h2   = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                          tmp2.x, tmp2.y, tmp2.z)],
                                                      mesh->serialMeshData[mesh->templNumb(
                                                          tmp1.x, tmp1.y, tmp1.z)]);
                                        v = (h1 - h_) * (v1 - v2) / h2 + v1;
                                    }
                                    tmp_V = tmp_V + v;
                                }
                            }
                            tmp1 = mesh->doStepY(tmp, 1);
                            points_edge.point2 =
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                            if ((mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                     .isOut == false))
                            {
                                if (isIntersection(bp.edge, points_edge, eps_p))
                                {
                                    h_ = distToEdgeY(bp.edge, mesh->serialMeshData[mesh->templNumb(
                                                                  tmp1.x, tmp1.y, tmp1.z)]);
                                    h1 = distP2PY(
                                        mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                        mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y,
                                                                             tmp1.z)]);
                                    v1 = V[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                                    v2 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                    v  = h1 * (v1 - v2) / h_ + v2;
                                    nm++;
                                    if (h_ < eps_h)
                                    {
                                        tmp2 = mesh->doStepY(tmp1, 1);
                                        v2   = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                        h2   = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                          tmp2.x, tmp2.y, tmp2.z)],
                                                      mesh->serialMeshData[mesh->templNumb(
                                                          tmp1.x, tmp1.y, tmp1.z)]);
                                        v = (h1 - h_) * (v1 - v2) / h2 + v1;
                                    }
                                    tmp_V = tmp_V + v;
                                }
                            }

                            tmp1 = mesh->doStepY(tmp, -1);
                            points_edge.point2 =
                                mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                            if ((mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                     .isOut == false))
                            {
                                if (isIntersection(bp.edge, points_edge, eps_p))
                                {
                                    h_ = distToEdgeY(bp.edge, mesh->serialMeshData[mesh->templNumb(
                                                                  tmp1.x, tmp1.y, tmp1.z)]);
                                    h1 = distP2PY(
                                        mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)],
                                        mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y,
                                                                             tmp1.z)]);
                                    v1 = V[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                                    v2 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                                    v  = h1 * (v1 - v2) / h_ + v2;
                                    nm++;
                                    if (h_ < eps_h)
                                    {
                                        tmp2 = mesh->doStepY(tmp1, -1);
                                        v2   = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                                        h2   = distP2PX(mesh->serialMeshData[mesh->templNumb(
                                                          tmp2.x, tmp2.y, tmp2.z)],
                                                      mesh->serialMeshData[mesh->templNumb(
                                                          tmp1.x, tmp1.y, tmp1.z)]);
                                        v = (h1 - h_) * (v1 - v2) / h2 + v1;
                                    }
                                    tmp_V = tmp_V + v;
                                }
                            }

                            if (nm != 0)
                            {
                                V[mesh->templNumb(tmp.x, tmp.y, tmp.z)] = tmp_V / (1.0 * nm);
                            }
                            else
                            {
                                strange_points.push_back(tmp);
                            }
                        }
                    }
                }
            }
        }
    }

    PointType                  v3;
    DGeo::Point<int>           tmp3;
    BoundaryPoint3d<PointType> bp_1, bp_2;
    for (int g = 0; g < 10; ++g)
    {
        for (int i = 0; i < strange_points.size(); ++i)
        {
            tmp_V = 0;
            nm    = 0;
            tmp   = strange_points[i];
            // 1
            tmp1 = mesh->doStepX(tmp, 1);
            tmp2 = mesh->doStepY(tmp1, 1);
            tmp3 = mesh->doStepY(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {

                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 2
            tmp1 = mesh->doStepX(tmp, 1);
            tmp2 = mesh->doStepY(tmp1, -1);
            tmp3 = mesh->doStepY(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 3
            tmp1 = mesh->doStepY(tmp, -1);
            tmp2 = mesh->doStepX(tmp1, -1);
            tmp3 = mesh->doStepX(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 4
            tmp1 = mesh->doStepX(tmp, -1);
            tmp2 = mesh->doStepY(tmp1, 1);
            tmp3 = mesh->doStepY(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 5
            tmp1 = mesh->doStepZ(tmp, 1);
            tmp2 = mesh->doStepY(tmp1, 1);
            tmp3 = mesh->doStepY(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 6
            tmp1 = mesh->doStepZ(tmp, 1);
            tmp2 = mesh->doStepY(tmp1, -1);
            tmp3 = mesh->doStepY(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 7
            tmp1 = mesh->doStepY(tmp, -1);
            tmp2 = mesh->doStepZ(tmp1, -1);
            tmp3 = mesh->doStepZ(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 8
            tmp1 = mesh->doStepZ(tmp, -1);
            tmp2 = mesh->doStepY(tmp1, 1);
            tmp3 = mesh->doStepY(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 9
            tmp1 = mesh->doStepX(tmp, 1);
            tmp2 = mesh->doStepZ(tmp1, 1);
            tmp3 = mesh->doStepZ(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 10
            tmp1 = mesh->doStepX(tmp, 1);
            tmp2 = mesh->doStepZ(tmp1, -1);
            tmp3 = mesh->doStepZ(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 11
            tmp1 = mesh->doStepZ(tmp, -1);
            tmp2 = mesh->doStepX(tmp1, -1);
            tmp3 = mesh->doStepX(tmp, -1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            // 12
            tmp1 = mesh->doStepX(tmp, -1);
            tmp2 = mesh->doStepZ(tmp1, 1);
            tmp3 = mesh->doStepZ(tmp, 1);
            if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y) || (tmp1.z != tmp.z)) &&
                ((tmp3.x != tmp.x) || (tmp3.y != tmp.y) || (tmp3.z != tmp.z)) &&
                ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y) || (tmp2.z != tmp1.z)))
            {
                v1    = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                v2    = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];
                v3    = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                tmp_V = tmp_V + v1 + v3 - v2;
                nm++;
            }
            V[mesh->templNumb(tmp.x, tmp.y, tmp.z)] = tmp_V / (1.0 * nm);
        }
    }

    strange_points.clear();
};

template <class PointType>
bool PoissonSolver3d<PointType>::isIntersection(DGeo::Edge<PointType> boundary_edge,
                                                DGeo::Edge<PointType> points_edge, double eps)
{
    if (polar == false)
    {
        return boundary_edge.IsIntersectionEdge(points_edge, eps);
    }
    else
    {
        DGeo::Edge<PointType> edge;
        Dmath::Polar2Cartesian(points_edge.point1.x, points_edge.point1.y, edge.point1.x,
                               edge.point1.y);
        Dmath::Polar2Cartesian(points_edge.point2.x, points_edge.point2.y, edge.point2.x,
                               edge.point2.y);
        return boundary_edge.IsIntersectionEdge(edge, eps);
    }
};

template <class PointType>
void PoissonSolver3d<PointType>::fieldCalculation_Simple_FDE_3d(
    std::vector<PointType>& V, const std::shared_ptr<MeshContainer3d<PointType>>& mesh, char type,
    std::vector<PointType>& E)
{

    double                     phi0, phi1, phi2, phi3;
    DGeo::Point<int>           tmp, tmp1, tmp2, tmp3;
    BoundaryPoint3d<PointType> bp;
    double                     h_, h1, h2, h3, delta, H;

    int flagE = 0;
    int m     = -1;
    for (int k = 0; k < mesh->templNumb.GetNz() - 1; k++)
    {
        tmp1.z = k;
        for (int i = 0; i < mesh->templNumb.GetNrow() - 1; i++)
        {
            tmp1.x = i;
            for (int j = 0; j < mesh->templNumb.GetNcol() - 1; j++)
            {
                tmp1.y = j;
                if (mesh->templNumb(i, j, k) != -1)
                {
                    ++m;

                    tmp.x = j;
                    tmp.y = i;
                    tmp.z = k;

                    flagE = 0;
                    tmp1  = doStep(mesh, type, tmp, -1);
                    tmp2  = doStep(mesh, type, tmp, 1);

                    phi1 = V[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)];
                    phi0 = V[mesh->templNumb(tmp.x, tmp.y, tmp.z)];
                    phi2 = V[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)];

                    h1 = distP2P(mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                 mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)], type);
                    h2 = distP2P(mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)],
                                 mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)], type);

                    if ((tmp.x == tmp1.x) && (tmp1.y == tmp.y) && (tmp1.z == tmp.z) ||
                        (mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)].isOut ==
                         true) &&
                            (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)].isOut ==
                             true) &&
                            (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)].isOut ==
                             false))
                    {
                        // 1
                        // tmp,tmp1 - 1, tmp2 - 2
                        bp   = findBoundaryPoint(tmp);
                        tmp3 = doStep(mesh, type, tmp2, 1);
                        phi3 = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                        h3 = distP2P(mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)],
                                     mesh->serialMeshData[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)],
                                     type);
                        H     = h2 + h3;
                        delta = h3 / h2;
                        E[m]  = -(-(2 + delta) * phi0 + (1 + delta) * (1 + delta) * phi2 / delta -
                                 phi3 / delta) /
                               H;
                    }
                    else if ((tmp.x == tmp2.x) && (tmp.y == tmp2.y) && (tmp.z == tmp2.z) ||
                             (mesh->serialMeshData[mesh->templNumb(tmp.x, tmp.y, tmp.z)].isOut ==
                              true) &&
                                 (mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)]
                                      .isOut == false) &&
                                 (mesh->serialMeshData[mesh->templNumb(tmp2.x, tmp2.y, tmp2.z)]
                                      .isOut == true))
                    {
                        // 5
                        // tmp,tmp2 - 5, tmp1 - 4
                        tmp3 = doStep(mesh, type, tmp1, -1);
                        phi3 = V[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)];
                        h3 = distP2P(mesh->serialMeshData[mesh->templNumb(tmp1.x, tmp1.y, tmp1.z)],
                                     mesh->serialMeshData[mesh->templNumb(tmp3.x, tmp3.y, tmp3.z)],
                                     type);
                        H     = h1 + h3;
                        delta = h1 / h3;
                        E[m]  = -(delta * phi3 - (1 + delta) * (1 + delta) * phi1 / delta +
                                 (2 + delta) * phi0 / delta) /
                               H;
                    }
                    else
                    {
                        // 3
                        // tmp1 - 2, tmp2 - 4
                        H     = h1 + h2;
                        delta = h2 / h1;
                        E[m] =
                            -(-delta * phi1 + (delta * delta - 1) * phi0 / delta + phi2 / delta) /
                            H;
                    }
                }
            }
        }
    }
};

template <class PointType>
double PoissonSolver3d<PointType>::distToEdge(DGeo::Edge<PointType>  edge,
                                              DGeo::Point<PointType> point, char type)
{
    if (type == 'x')
    {
        return distToEdgeX(edge, point);
    }
    else
    {
        return distToEdgeY(edge, point);
    }
};

template <class PointType>
DGeo::Point<int>
PoissonSolver3d<PointType>::doStep(const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                                   char type, DGeo::Point<int> point, int step)
{
    if (type == 'x')
    {
        return mesh->doStepX(point, step);
    }
    else if (type == 'y')
    {
        return mesh->doStepY(point, step);
    }
    else
    {
        return mesh->doStepZ(point, step);
    }
};
