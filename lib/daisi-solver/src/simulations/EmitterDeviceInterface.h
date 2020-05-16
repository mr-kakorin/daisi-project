#ifndef EMITTERDEVICEINTERFACE_H
#define EMITTERDEVICEINTERFACE_H
#include "BoundaryContainer2d.h"
#include "ParticleSource.h"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include <algorithm>
#include <vector>
template <class PointType, class ParticlesDataType>
class EmitterDeviceInterface
{
    friend class boost::serialization::access;

  protected:
    int              emissionType;
    int              distributionEnergyType; // ��� ��������������� �������������
    int              generateType;      // ��� ������������� (���������/�����������)
    int              nParticlesEnergy;  // ����� ������������ ������ � ��������� ��������
    double           energyAverage;     // ������� ������� ������������ ������
    double           energySpread;      // �������������� ������c ������������ ������
    int              nParticlesEmitter; //����� ������������ ������ � ������ ���������� �� ��������
    double           totalCurrent;
    std::vector<int> boundaryList;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& emissionType;
        ar& distributionEnergyType;
        ar& generateType;
        ar& nParticlesEnergy;
        ar& energyAverage;
        ar& energySpread;
        ar& nParticlesEmitter;
        ar& boundaryList;
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& emissionType;
        ar& distributionEnergyType;
        ar& generateType;
        ar& nParticlesEnergy;
        ar& energyAverage;
        ar& energySpread;
        ar& nParticlesEmitter;
        ar& boundaryList;
    };

  public:
    std::vector<unsigned int> newIndexes;

    virtual std::vector<std::vector<float>> GetCurrentDensityDistribution() = 0;

    std::vector<int> GetBoundariesList()
    {
        return boundaryList;
    };
    void SetBoundariesList(std::vector<int>                            in,
                           std::vector<BoundaryContainer2d<PointType>> boundaryIn)
    {
        boundaryList = in;
        InitEmissionBoundary(boundaryIn);
    };
    virtual void InitEmissionBoundary(std::vector<BoundaryContainer2d<PointType>> boundaryIn) = 0;
    EmitterDeviceInterface()
    {
        emissionType           = 0;
        distributionEnergyType = 0;
        generateType           = 0;
        nParticlesEnergy       = 1;
        nParticlesEmitter      = 1;
        energyAverage          = 0;
        energySpread           = 0;
    };
    void SetEmissionType(int type)
    {
        emissionType = type;
    };
    int GetEmissionType()
    {
        return emissionType;
    };
    virtual std::vector<int> GetParticlesNumber()                  = 0;
    virtual void             SetParticlesNumber(std::vector<int>)  = 0;
    virtual void SetDistributionParameters(std::vector<double> in) = 0;
    virtual std::vector<double> GetDistributionParameters()        = 0;
    virtual double              GetEmissionCurrent()               = 0;
    virtual void SetEmissionCurrent(double current)                = 0;
    virtual void GenerateParticles(const std::vector<unsigned int>& EmptyPlaces,
                                   ParticlesDataType& particlesData, PointType restMass,
                                   int flagClear) = 0;
    virtual int GetNumbersOfParticlesGeneration() = 0;
};
#endif