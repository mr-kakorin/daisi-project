#ifndef PARTICLESMOVERBASE_H
#define PARTICLESMOVERBASE_H
#include <boost/archive/binary_iarchive.hpp>
#include <Constants.h>

template <class PointType>
class Particles2dpolar;
template <class PointType>
class Particles3d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class Particles3dcil;
template <class PointType>
class Particles2d;
template <class PointType>
class particlesFields;
template <class PointType>

class ParticlesMover
{
  private:
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    //	int moverType;

    //	PointType timeStep;

    std::vector<std::vector<double>> timeSteps;

    std::vector<PointType> tmp1;
    std::vector<PointType> tmp2;
    std::vector<PointType> DT;

  public:
    std::vector<double> params;
    template <class particlesdataType>
    void stepEstimate(const std::shared_ptr<particlesdataType>& particlesData, int nlow, int i1,
                      int i2, int thread, int solverType);

    void InitParallel(int numThreads);
    int flagInit;

    void init(std::vector<PointType> alpha, int solverType);
    ParticlesMover();
    PointType GetTimeStep(int i, int thread);
    void SetTimeSteps(std::vector<double> in);
    void SetParameters(std::vector<double> in);
    void                DoubleTimeSteps();
    void                HalveTimeSteps();
    std::vector<double> GetParameters();
    void updatePositions(const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                         const std::shared_ptr<Particles3dcil<PointType>>& particlesDataTmp,
                         int nlow, int i1, int i2, int thread);
    void updateMomentums(const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                         particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow,
                         int i1, int i2, int thread);
    void updatePositions(const std::shared_ptr<Particles3d<PointType>>& particlesData,
                         const std::shared_ptr<Particles3d<PointType>>& particlesDataTmp, int nlow,
                         int i1, int i2, int thread);
    void updateMomentums(const std::shared_ptr<Particles3d<PointType>>& particlesData,
                         particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow,
                         int i1, int i2, int thread);
    void updatePositions(const std::shared_ptr<Particles2d<PointType>>& particlesData,
                         const std::shared_ptr<Particles2d<PointType>>& particlesDataTmp, int nlow,
                         int i1, int i2, int thread);
    void updateMomentums(const std::shared_ptr<Particles2d<PointType>>& particlesData,
                         particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow,
                         int i1, int i2, int thread);
    void updatePositions(const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
                         const std::shared_ptr<Particles2dpolar<PointType>>& particlesDataTmp,
                         int nlow, int i1, int i2, int thread);
    void updateMomentums(const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
                         particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow,
                         int i1, int i2, int thread);
};
#endif