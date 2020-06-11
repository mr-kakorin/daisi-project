#ifndef MESHGENERATOR2D_H
#define MESHGENERATOR2D_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <thread>

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
}

namespace Dmath
{
class imat;
}

template <class PointType>
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;
template <class PointType>
class MeshContainer3d;
template <class PointType>
class MeshContainer2d;
template <class PointType>
class ParticleGridInterface;

class MeshParams
{
  public:
    std::vector<double> leftBorder[3];

    std::vector<double> rightBorder[3];

    std::vector<int> type[3];

    std::vector<std::vector<double>> params[3];

    MeshParams(std::string filename, double& eps, std::string& errorMsg);
    double GetStep(double point, int flagX) const noexcept;
    int GetNumberOfSteps(int flagX, double min, double max, double& hMin) const;
};

template <class PointType>
class MeshGenerator
{
    friend class boost::serialization::access;

  private:
    bool                flagBoundaryInit;
    bool                flagMeshInit;
    std::vector<int>    boundaryList;
    std::vector<double> meshParam;
    int                 t1;
    int                 t2;

    void RayTrace(double epsilon, const MeshParams& meshParams, DGeo::Point<PointType> StartPoint,
                  Dmath::imat& flagMartix, PointType h1, PointType h2, int yNum,
                  const std::shared_ptr<BoundaryContainer2d<PointType>> boundary);
    void RayTracePolar(double epsilon, const MeshParams& meshParams,
                       DGeo::Point<PointType> StartPoint, Dmath::imat& flagMartix, PointType h1,
                       PointType h2, int yNum,
                       const std::shared_ptr<BoundaryContainer2d<PointType>> boundary);

    void MeshAssembly(const MeshParams& meshParams, DGeo::Point<PointType>,
                      std::vector<DGeo::Point<PointType>>*, Dmath::imat&, Dmath::imat&, PointType,
                      PointType, int, int&,
                      const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                      std::vector<int>&                                     linearToTemplNumb);
    void MeshAssemblyPolar(const MeshParams& meshParams, DGeo::Point<PointType>,
                           std::vector<DGeo::Point<PointType>>*, Dmath::imat&, Dmath::imat&,
                           PointType, PointType, int, int&,
                           const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                           std::vector<int>&                                     linearToTemplNumb);

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    std::vector<double> GetMeshParam();
    void SetMeshParam(std::vector<double> in);
    std::vector<int> GetBoundaryList();
    void SetBoundaryList(std::vector<int>                            in,
                         std::vector<BoundaryContainer2d<PointType>> boundaryIn);
    void Convert2ParticleGridInterface(ParticleGridInterface<PointType>*);

    bool IsInitBoundary();
    bool IsMeshGen();
    MeshGenerator(std::string);
    MeshGenerator();
    void WriteMesh2VTK(std::string InputFileName){};
    void WriteBoundary2VTK(std::string InputFileName){};
    void MeshGenerate(double epsilon, const MeshParams& meshParams, double&,
                      const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                      const std::shared_ptr<MeshContainer2d<PointType>>&    mesh,
                      std::string&                                          errorMsg);
    void MeshGenerate(std::string                                           meshParam, double&,
                      const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                      const std::shared_ptr<MeshContainer2d<PointType>>& mesh, int flagVTK,
                      std::string& errorMsg);

    void MeshGenerate(std::string meshParamFile, double& progress,
                      const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                      const std::shared_ptr<MeshContainer3d<PointType>>& mesh, int flagVTK,
                      std::string& errorMsg);
    void MeshGenerate(std::string meshParamFile, double& progress,
                      const std::shared_ptr<BoundaryContainer3d<PointType>> boundary,
                      const std::shared_ptr<MeshContainer3d<PointType>>& mesh, int flagVTK,
                      std::string& errorMsg);

    void MeshGeneratePolar(std::string                                           meshParam, double&,
                           const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                           const std::shared_ptr<MeshContainer2d<PointType>>& mesh, int flagVTK,
                           std::string& errorMsg);
    void MeshGeneratePolar(double epsilon, const MeshParams& meshParams, double&,
                           const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
                           const std::shared_ptr<MeshContainer2d<PointType>>&    mesh,
                           std::string&                                          errorMsg);
};

#endif