#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H
#include "PoissonSolverBase.h"
#include <boost/archive/binary_iarchive.hpp>
#include <cmath>

typedef double (*coeff)(double, double, double, double, double, double);
typedef double (*coeff3d)(double, double, double, double, double, double, double, double, double);

template <class PointType>
class MeshContainer2d;
template <class PointType>
class MeshContainer3d;
template <class PointType>
class GridData2d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
class BoundaryConditions;

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
}

template <class PointType>
class PoissonSolver
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void                load(Archive& ar, const unsigned int);
    std::vector<double> param;

  public:
    bool                                  polar;
    int                                   ChargeSpace;
    int                                   RecalculateParameter;
    double                                w_charge;
    double                                eps_tolerance_charge;
    std::vector<int>                      kZeros;
    std::vector<int>                      solverFlags;
    std::vector<BoundaryPoint<PointType>> boundaryPoints;
    std::vector<IBPoint<PointType>>       ibPoints;
    std::vector<int>                      nearBoundaryPoints_n;
    std::vector<int>                      nearBoundaryPoints_all;
    std::vector<int>                      rights;
    std::vector<int>                      lefts;
    std::vector<int>                      middles;
    std::vector<int>                      ups;
    std::vector<int>                      downs;
    std::vector<int>                      deeps;
    std::vector<int>                      outwards;
    std::vector<double>                   c_rights;
    std::vector<double>                   c_lefts;
    std::vector<double>                   c_ups;
    std::vector<double>                   c_downs;
    std::vector<double>                   c_middles;
    std::vector<double>                   c_deeps;
    std::vector<double>                   c_outwards;
    std::vector<double>                   vect;
    std::vector<PointType>                xCharge;
    std::vector<PointType>                linVect;
    std::vector<PointType>                x;

    std::vector<BoundaryPointIntersection<PointType>> boundaryPointsIntersectionDefault;
    std::vector<BoundaryPointIntersection<PointType>> boundaryPointsIntersectionCondition;

    std::vector<NonZeroLinVect<PointType>> nonZeroLinVect;

    double    w;
    double    eps_tolerance;
    PointType eps_p; // for equals of 2 points
    PointType eps_h; // for cut cells

    int nSlicePoints;
    int nSliceBPoints;

    PoissonSolver();

    void                SetParameters(const std::vector<double>& param);
    std::vector<double> GetParameters();

    int                                               systemSize;
    int                                               boundCondNum;
    std::shared_ptr<MeshContainer2d<PointType>> mesh_copy;

    void searchBoundaryPointsIntersection(
        const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
        const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, double& progress);

    DGeo::Point<int> closestVertex(const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                                   DGeo::Point<PointType>                             point);

    double dist(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distToEdgeX(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgeY(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgeZ(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgePolarPhi(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdgePolarR(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point);
    double distToEdge(DGeo::Edge<PointType> edge, DGeo::Point<PointType> point, char type);
    BoundaryPoint<PointType> findBoundaryPoint(int pnumber);
    DGeo::Point<int> getImagineMeshPoint(const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                                         DGeo::Point<int>                                   mp);
    void             addNearBoundaryPoint_n(int n);
    void             addNearBoundaryPoint_all(int n);
    void             CreateSolutionSpace(const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                                         std::vector<int>&                                  nonZeros);

    std::vector<int> getNearBoundarypoints_n();
    std::vector<int> getNearBoundarypoints_all();

    std::vector<int> nonZerosV;

    std::vector<int>       indexes;
    std::vector<int>       templ;
    std::vector<int>       n2k;
    std::vector<PointType> Vtmp;

    void addBoundaryPoint(BoundaryPoint<PointType> bpoint, DGeo::Point<PointType> ipoint);
    void addBoundaryPoint(BoundaryPoint<PointType> bpoint);
    void noteBoundaryOnMesh(
        const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
        const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, double& progress);
    void createMatrix(const std::shared_ptr<MeshContainer2d<PointType>>& mesh, coeff left,
                      coeff right, coeff up, coeff down, coeff middle);
    void createMatrix3d(const std::shared_ptr<MeshContainer3d<PointType>>& mesh, coeff3d left,
                        coeff3d right, coeff3d up, coeff3d down, coeff3d outward, coeff3d deep,
                        coeff3d middle);
    void InitSolver(const std::shared_ptr<GridData2d<PointType>>&                       gridData,
                    const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions,
                    double&                                          progress);
    void InitSolver(const std::shared_ptr<GridData2daxs<PointType>>&                    gridData,
                    const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions,
                    double&                                          progress);
    void InitSolver(const std::shared_ptr<GridData2dpolar<PointType>>&                  gridData,
                    const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions,
                    double&                                          progress);
    void InitSolver(const std::shared_ptr<GridData3d<PointType>>&                       gridData,
                    const std::shared_ptr<MeshContainer3d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions,
                    double&                                          progress);
    void InitSolver(const std::shared_ptr<GridData3d<PointType>>&                       gridData,
                    const std::shared_ptr<MeshContainer3d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions,
                    double&                                          progress){

    };

    void FieldSimulate(const std::shared_ptr<GridData2d<PointType>>&      gridData,
                       const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                       const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                       double t, double& progress);
    void FieldSimulate(const std::shared_ptr<GridData2daxs<PointType>>&   gridData,
                       const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                       const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                       double t, double& progress);
    void FieldSimulate(const std::shared_ptr<GridData2dpolar<PointType>>& gridData,
                       const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                       const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                       double t, double& progress);
    void FieldSimulate(const std::shared_ptr<GridData3d<PointType>>&      gridData,
                       const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                       const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                       double t, double& progress);

    void FieldSimulateCharge(const std::shared_ptr<GridData2d<PointType>>&      gridData,
                             const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                             const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                             double t, std::vector<int>& nonZeros, int step);
    void FieldSimulateCharge(const std::shared_ptr<GridData2daxs<PointType>>&   gridData,
                             const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                             const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                             double t, std::vector<int>& nonZeros, int step);
    void FieldSimulateCharge(const std::shared_ptr<GridData2dpolar<PointType>>& gridData,
                             const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                             const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                             double t, std::vector<int>& nonZeros, int step){

    };
    void FieldSimulateCharge(const std::shared_ptr<GridData3d<PointType>>&      gridData,
                             const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                             const std::shared_ptr<BoundaryConditions>&   boundaryConditions,
                             double t, std::vector<int>& nonZeros, int step){

    };

    void fieldCalculation_Simple_FDE(std::vector<PointType>& rho, std::vector<PointType>& V,
                                     const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                                     char type, std::vector<PointType>& E);
    void fieldCalculation_Simple_FDE(std::vector<PointType>& rho, std::vector<PointType>& V,
                                     const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                                     char type, std::vector<PointType>& E);
    DGeo::Point<int> doStep(const std::shared_ptr<MeshContainer2d<PointType>>& mesh, char type,
                            DGeo::Point<int> point, int step);
    DGeo::Point<int> doStep(const std::shared_ptr<MeshContainer3d<PointType>>& mesh, char type,
                            DGeo::Point<int> point, int step);
    void             searchBoundaryPointsIntersectionPolar(
                    const std::shared_ptr<MeshContainer2d<PointType>>&                  mesh,
                    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions, double& progress);
    void solve(std::vector<PointType>& rho, std::vector<PointType>& V,
               const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
               const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t,
               double& progress);
    void solveCharge(std::vector<PointType>& rho, std::vector<PointType>& V,
                     const std::shared_ptr<MeshContainer2d<PointType>>& mesh,
                     const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t,
                     std::vector<int>& nonZeros);

    void   solve_3d(std::vector<PointType>& rho, std::vector<PointType>& V,
                    const std::shared_ptr<MeshContainer3d<PointType>>& mesh,
                    const std::shared_ptr<BoundaryConditions>& boundaryConditions, double t);
    double distP2PX(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distP2PY(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    double distP2PZ(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2);
    bool   isIntersection(DGeo::Edge<PointType> boundary_edge, DGeo::Edge<PointType> points_edge,
                          double eps);
    double distP2P(DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2, char type);
};

#endif
