#ifndef EMCS_H
#define EMCS_H

#include "memory"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
template <class PointType>
struct CurvePoints_t;
};

template <class PointType>
class EmitterDevice2daxs;
template <class PointType>
class EmitterDevice2d;
template <class PointType>
class EmitterDevice3d;
template <class PointType>
class ParticleSource2d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class GridData2d;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
template <class PointType>
class GridData3dpolar;
template <class PointType>
class ParticleGridInterface;
template <class PointType>
class NearCathodeVolume;
template <class PointType>

class EmissionCurrentSolverBase
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    int                                                    algorithm;
    std::vector<double>                                    Lem;
    std::vector<double>                                    Hem;
    std::vector<int>                                       flowsNumbers;
    std::vector<std::vector<unsigned int>>                 cellNumbers;
    std::vector<unsigned int>                              emissionCells;
    std::vector<std::vector<int>>                          istarts;
    std::vector<std::vector<DGeo::Point<PointType>>>       points1;
    std::vector<std::vector<DGeo::Point<PointType>>>       points2;
    std::vector<DGeo::CurvePoints_t<PointType>>            gradients;
    std::vector<std::vector<DGeo::Edge<PointType>>>        emittingCutEdge;
    std::vector<std::vector<double>>                       CathodeFields;
    std::vector<std::vector<NearCathodeVolume<PointType>>> nearCathodeVolumes;
    std::vector<std::vector<double>>                       E0;
    std::vector<std::vector<double>>                       K;

    template <class gridDataType>
    void CalculateCathodeFields(const std::shared_ptr<ParticleSource2d<PointType>>& source,
                                const std::shared_ptr<gridDataType>& gridData, int flowNumber);

    void CalculateCathodeFields(const std::shared_ptr<ParticleSource2d<PointType>>& source,
                                const std::shared_ptr<GridData3d<PointType>>&       gridData,
                                int                                                 flowNumber){

    };

    void SetValueOnSource(const std::shared_ptr<ParticleSource2d<PointType>>& source,
                          std::vector<double> value, int flowNumber, int flag);

    void reset();

    template <class gridDataType>
    void init(int i, const std::shared_ptr<gridDataType>& gridData, int flagEm,
              const std::shared_ptr<ParticleSource2d<PointType>>&      source,
              const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
              int                                                      flag);

    void CalculateCathodeFields(std::vector<std::shared_ptr<ParticleSource2d<PointType>>>& source,
                                const std::shared_ptr<GridData3d<PointType>>&              gridData,
                                int flowNumber);
    void init(int i, const std::shared_ptr<GridData3d<PointType>>& gridData, int flagEm,
              std::vector<std::shared_ptr<ParticleSource2d<PointType>>> source,
              const std::shared_ptr<ParticleGridInterface<PointType>>&  particleGridInterface,
              int                                                       flag);
    int  GetEmSize();
    void SetParameters(std::vector<std::vector<double>> in);
    std::vector<std::vector<double>> GetParameters();
    void addFlow(int flowNumber);
    EmissionCurrentSolverBase();
};
template <class PointType>
class EmissionCurrentSolverPTI : public EmissionCurrentSolverBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<EmissionCurrentSolverBase<PointType>>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<EmissionCurrentSolverBase<PointType>>(*this);
    }
    PointType polinomCalc(const std::vector<PointType>& polinom, PointType x);

  public:
    /*std::vector<PointType> GetPolinom()
    {
            return polinom;
    };*/
    /*std::vector<unsigned int> GetEmissionCells()
    {
            return emissionCells;
    };
    void ChargeConserving(const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge);





    std::vector<double> CalcCathodeFieldsPTI( const std::shared_ptr<ParticleSource2d<PointType>>&
    source, const std::shared_ptr<GridData2daxs<PointType>>& gridData, int nE, int flowNumber);
    std::vector<double> CalcCathodeFieldsPTI( const std::shared_ptr<ParticleSource2d<PointType>>&
    source, const std::shared_ptr<GridData2d<PointType>>& gridData, int nE, int flowNumber);
    std::vector<double> CalcCathodeFieldsPTI( const std::shared_ptr<ParticleSource2d<PointType>>&
    source, const std::shared_ptr<GridData2dpolar<PointType>>& gridData, int nE, int flowNumber);

    void SetEmissionCurrent(const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    PointType current); void SetEmissionCurrent(const
    const std::shared_ptr<EmitterDevice2d<PointType>>& emitter, PointType current); void
    SetEmissionCurrent(const std::shared_ptr<EmitterDevice3d<PointType>>& emitter, PointType
    current);

    void UpdateEmissionCurrent(int& flagE, const std::shared_ptr<EmitterDevice2d<PointType>>&
    emitter, const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface, const
    const std::shared_ptr<GridData2d<PointType>>& gridData, PointType timeStep, int flowNumber, int
    stepNumber, double mass, double charge); void UpdateEmissionCurrent(int& flagE,
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge); void UpdateEmissionCurrent(int& flagE, const
    const std::shared_ptr<EmitterDevice3d<PointType>>& emitter, const
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface, const
    const std::shared_ptr<GridData3d<PointType>>& gridData, PointType timeStep, int flowNumber, int
    stepNumber, double mass, double charge); void UpdateEmissionCurrent(int& flagE,
    const std::shared_ptr<EmitterDevice2d<PointType>>& emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge);






            template <class gridDataType, class emittersType>
    void UpdateEmissionCurrent(int& flagE, arma::mat Icoef, emittersType emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<gridDataType>& gridData, int flagStopUpdateCurrent, std::vector<double>
    ChargeSign);

    void UpdateEmissionCurrent(int& flagE, arma::mat Icoef,
    std::vector<std::shared_ptr<EmitterDevice3d<PointType>>&> emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<GridData3d<PointType>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);


    void VirtualDiode(int& flagE, const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge); void NewtonIterative(const
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter, const
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface, const
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int
    stepNumber, double mass, double charge);

    void ChangePolinom(const std::shared_ptr<EmitterDevice2d<PointType>>& emitter,
    std::vector<double> dJ,int p1,int p2)
    {

    };
    void ChangePolinom(const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    std::vector<double> dJ, int p1, int p2);

    double ChangeCurrent(const std::shared_ptr<EmitterDevice2d<PointType>>& emitter, int i, double
    per)
    {
            return 0;
    };
    double ChangeCurrent(const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter, int i,
    double per);*/
};
template <class PointType>
class EmissionCurrentSolverPIC : public EmissionCurrentSolverBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<EmissionCurrentSolverBase<PointType>>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<EmissionCurrentSolverBase<PointType>>(*this);
    }

  public:
    void SetEmissionCurrent(const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
                            PointType                                             current);
    void SetEmissionCurrent(const std::shared_ptr<EmitterDevice2d<PointType>>& emitter,
                            PointType                                          current);
    void SetEmissionCurrent(const std::shared_ptr<EmitterDevice3d<PointType>>& emitter,
                            PointType                                          current);
    void UpdateEmissionCurrent(
        const std::shared_ptr<EmitterDevice2d<PointType>>&       emitter,
        const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
        const std::shared_ptr<GridData2d<PointType>>& gridData, PointType timeStep, int flowNumber,
        int stepNumber, double mass, double charge);
    void UpdateEmissionCurrent(
        const std::shared_ptr<EmitterDevice2d<PointType>>&       emitter,
        const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
        const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType timeStep,
        int flowNumber, int stepNumber, double mass, double charge);
    void UpdateEmissionCurrent(
        const std::shared_ptr<EmitterDevice2daxs<PointType>>&    emitter,
        const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
        const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep,
        int flowNumber, int stepNumber, double mass, double charge);
    void UpdateEmissionCurrent(
        const std::shared_ptr<EmitterDevice3d<PointType>>&       emitter,
        const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
        const std::shared_ptr<GridData3d<PointType>>& gridData, PointType timeStep, int flowNumber,
        int stepNumber, double mass, double charge);

    template <class gridDataType, class emittersType>
    void UpdateEmissionCurrent(
        const std::shared_ptr<emittersType>&                     emitter,
        const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
        const std::shared_ptr<gridDataType>& gridData, double timeStep, int flowNumber,
        int stepNumber, double mass, double charge, int emissionType);

    void
    VirtualDiode(const std::shared_ptr<EmitterDevice2daxs<PointType>>&    emitter,
                 const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
                 const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep,
                 int flowNumber, int stepNumber, double mass, double charge);
	void
	VirtualDiode1(const std::shared_ptr<EmitterDevice2daxs<PointType>>&    emitter,
			const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
			const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep,
			int flowNumber, int stepNumber, double mass, double charge);
    void
    ChargeConserving(const std::shared_ptr<EmitterDevice2daxs<PointType>>&    emitter,
                     const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
                     const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep,
                     int flowNumber, int stepNumber, double mass, double charge);
    void Poisson(const std::shared_ptr<EmitterDevice2daxs<PointType>>&    emitter,
                 const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface,
                 const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep,
                 int flowNumber, int stepNumber, double mass, double charge);

    EmissionCurrentSolverPIC(){};
};
#endif