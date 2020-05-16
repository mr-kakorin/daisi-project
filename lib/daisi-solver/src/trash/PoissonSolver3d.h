#ifndef POISSONSOLVER3D_H
#define POISSONSOLVER3D_H
#include "PoissonSolver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <cmath>

typedef double (*coeff3d)(double, double, double, double, double, double, double, double, double);

template <class PointType>
class MeshContainer3d;
class BoundaryConditions;

template <class PointType>
class GridData2d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;

template <class PointType>
class BoundaryContainer3d;
template <class PointType>
class BoundaryContainer2d;

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
};

template <class PointType>
class BoundaryPoint3d
{
  public:
    int                    pnum;  // number of mesh point
    char                   type;  // type of boundary condition
    double                 value; // value of boundary condition
    DGeo::Edge<PointType>  edge;
    int                    bcnum; // number of boundary condition;
    DGeo::Point<PointType> ipoint;
    DGeo::Point<int>       coord;
};

template <class PointType>
class PoissonSolver3d
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& w;
        ar& eps_tolerance;
        ar& w_charge;
        ar& eps_tolerance_charge;
        ar& eps_p;
        ar& eps_h;
        ar& RecalculateParameter;
        ar& solverFlags;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& w;
        ar& eps_tolerance;
        ar& w_charge;
        ar& eps_tolerance_charge;
        ar& eps_p;
        ar& eps_h;
        ar& RecalculateParameter;
        ar& solverFlags;
        solverFlags.resize(4);
    }

  public:
    std::vector<int> solverFlags;

    int RecalculateParameter;

    bool polar;

    std::vector<std::vector<std::vector<BoundaryPoint3d<PointType>>>> boundaryPoints;
    std::vector<int>                                                  nearBoundaryPoints_n;
    std::vector<int>                                                  nearBoundaryPoints_all;
    std::vector<int>                                                  rights;
    std::vector<int>                                                  lefts;
    std::vector<int>                                                  middles;
    std::vector<int>                                                  ups;
    std::vector<int>                                                  downs;
    std::vector<int>                                                  deeps;
    std::vector<int>                                                  outwards;
    std::vector<double>                                               c_rights;
    std::vector<double>                                               c_lefts;
    std::vector<double>                                               c_ups;
    std::vector<double>                                               c_downs;
    std::vector<double>                                               c_middles;
    std::vector<double>                                               c_deeps;
    std::vector<double>                                               c_outwards;

    std::vector<std::vector<PointType>> U;

    std::vector<BoundaryPointIntersection<PointType>> boundaryPointsIntersectionDefault;
    std::vector<BoundaryPointIntersection<PointType>> boundaryPointsIntersectionCondition;

    std::vector<NonZeroLinVect<PointType>> nonZeroLinVect;

    double    w;
    double    eps_tolerance;
    PointType eps_p; // for equals of 2 points
    PointType eps_h; // for cut cells

    int    ChargeSpace;
    double w_charge;
    double eps_tolerance_charge;

    int              nSlicePoints;
    std::vector<int> nSliceBPoints;

    PoissonSolver3d()
    {
        w                    = 1.9;
        eps_tolerance        = 0.0001;
        w_charge             = 1.7;
        eps_tolerance_charge = 0.01;
        eps_p                = 0.00001;
        eps_h                = 0.00001;
        ChargeSpace          = 5;
        RecalculateParameter = 1;
        solverFlags.resize(4);
        solverFlags[0] = 1; // System solver
        solverFlags[1] = 1; // time depending flag
        solverFlags[2] = 1; // space-charge flag
        solverFlags[3] = 1; // fast initialization
    };

    void SetParameters(const std::vector<double>& param);
    std::vector<double> GetParameters();

    double* linVect; // right hand side
    double* x;       // vector of solution
    int     systemSize;
    int     boundCondNum;

    void searchBoundaryPointsIntersection_3d(
        const std::shared_ptr<MeshContainer3d<PointType>>&           mesh,
        std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaries,
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, double z);

    DGeo::Point<int> closestVertex_3d(const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                                      DGeo::Point<PointType>                             point);

    double dist(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distToEdgeX(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgeY(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgeZ(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdge(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point, char type);
    BoundaryPoint3d<PointType> findBoundaryPoint(DGeo::Point<int> bpoint);
    void addNearBoundaryPoint_n(int n);
    void addNearBoundaryPoint_all(int n);

    std::vector<int> getNearBoundarypoints_n();
    std::vector<int> getNearBoundarypoints_all();

    void addBoundaryPoint(BoundaryPoint3d<PointType> bpoint, DGeo::Point<PointType> ipoint);
    void addBoundaryPoint(BoundaryPoint3d<PointType> bpoint);

    void
    noteBoundaryOnMesh_3d(const std::shared_ptr<MeshContainer3d<PointType>>&            mesh,
                          std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                          const std::shared_ptr<BoundaryConditions>& boundaryConditions);
    void createMatrix3d(const std::shared_ptr<MeshContainer3d<PointType>>& mesh, coeff3d left,
                        coeff3d right, coeff3d up, coeff3d down, coeff3d outward, coeff3d deep,
                        coeff3d middle);
    void InitSolver(const std::shared_ptr<GridData3d<PointType>>&                 gridData,
                    const std::shared_ptr<MeshContainer3d<PointType>>&            mesh,
                    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions);
    void InitSolver(const std::shared_ptr<GridData3d<PointType>>&                 gridData,
                    const std::shared_ptr<MeshContainer3d<PointType>>&            mesh,
                    std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions){

    };

    void FieldSimulate(const std::shared_ptr<GridData3d<PointType>>&      gridData,
                       const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                       const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t);
    void FieldSimulateCharge(const std::shared_ptr<GridData3d<PointType>>&      gridData,
                             const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                             const std::shared_ptr<BoundaryConditions>&         boundaryConditions,
                             double t, std::vector<int>& nonZeros, int step){

    };

    void fieldCalculation_Simple_FDE_3d(std::vector<PointType>&                            V,
                                        const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                                        char type, std::vector<PointType>& E);
    DGeo::Point<int> doStep(const std::shared_ptr<MeshContainer3d<PointType>>& mesh, char type,
                            DGeo::Point<int> point, int step);
    void solve_one_rotating(std::vector<PointType>& rho, std::vector<PointType>& V,
                            const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                            const std::shared_ptr<BoundaryConditions>& boundaryConditions, int i,
                            double z_begin, double z_end);
    void solve_charge(std::vector<PointType>& rho, std::vector<PointType>& V,
                      const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                      const std::shared_ptr<BoundaryConditions>&         boundaryConditions);
    void solve_offset(std::vector<PointType>&                            V,
                      const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                      const std::shared_ptr<BoundaryConditions>&         boundaryConditions);
    double distP2PX(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distP2PY(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distP2PZ(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    bool isIntersection(DGeo::Edge<PointType> boundary_edge, DGeo::Edge<PointType> points_edge,
                        double eps);
    double distP2P(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2, char type);
    void extrapolate_potential(std::vector<PointType>&                            V,
                               const std::shared_ptr<MeshContainer3d<PointType>>& mesh);
    void sum_solutions(std::vector<std::vector<PointType>>& U, std::vector<PointType>& V,
                       const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t);
};

#endif
