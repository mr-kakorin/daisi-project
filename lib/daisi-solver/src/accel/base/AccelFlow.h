#ifndef AccelFlow_H
#define AccelFlow_H
//#include "ModelTemplate.h"
#include <armadillo>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <unordered_map>

class myunsorted_map;
class LinacDynamicsLong
{
  public:
    int    NparticlesLong;
    double longQ;

    std::vector<double> p;
    std::vector<double> z;
    std::vector<int>    cellNumber;
    std::vector<double> Ez;
    void clearandResize(int nParticles);
    std::vector<void*> GetData();
    int checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                          const std::vector<double>& LL, const std::vector<double>& MinumalRadii,
                          double lambda);
};
class SynchrotronDevice;
class LinacDynamicsAccel
{
  public:
    int    Nparticles;
    double longQ;
    void   Init();
    double CalculateTransmission(double Xaperture, double Yaperture, double ChannelRel = 1.0);
    void removeParticle(std::vector<unsigned int>& indexes);
    void searchBadParticle(std::vector<unsigned int>& Indexes);
    void removeParticle();

    std::vector<double> z;
    std::vector<double> dz;
    std::vector<double> x;
    std::vector<double> dx;
    std::vector<double> y;
    std::vector<double> dy;

    std::vector<double> zStart;
    std::vector<double> dzStart;
    std::vector<double> xStart;
    std::vector<double> dxStart;
    std::vector<double> yStart;
    std::vector<double> dyStart;
    std::vector<int>    isInPerture;
    std::vector<int>    cellNumber;
    // std::vector<double> Ez;
    void clearandResize(int nParticles);
    // std::vector<void*> GetData();
    // int checkBadParticles(int& cell, std::vector<std::vector<float>>&  OutParams, const
    // std::vector<double>&  LL, const std::vector<double>&  MinumalRadii, double lambda);
    std::vector<void*> GetData();
    ~LinacDynamicsAccel(){};
};

class SynchrotronDynamics
{
  public:
    int    Nparticles;
    double longQ;
    void   Init();
    double CalculateTransmission(double Xaperture, double Yaperture, double ChannelRel = 1.0);
    void removeParticle();

    std::vector<arma::vec> position;
    std::vector<arma::vec> positionStart;

    std::vector<int> isInPerture;
    void clearandResize(int nParticles);
    ~SynchrotronDynamics(){};
};

class LinacDynamicsTransv
{
  public:
    int                 NparticlesTransv;
    double              transvQ;
    std::vector<double> p;
    std::vector<double> z;
    std::vector<double> S11x;
    std::vector<double> S12x;
    std::vector<double> S22x;
    std::vector<double> S11y;
    std::vector<double> S12y;
    std::vector<double> S22y;
    std::vector<int>    cellNumber;
    std::vector<double> Ex;
    std::vector<double> Ey;
    std::vector<double> Ez;
    std::vector<double> EzCol;

    void clearandResize(int nParticles);
    std::vector<void*> GetData();
    int checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                          const std::vector<double>& LL, const std::vector<double>& MinumalRadii,
                          double lambda, double relChann);
};

class AccelFlowBase
{
  public:
    AccelFlowBase(int ProblemType);
    AccelFlowBase(){};
    void translateParameters();
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    myunsorted_map* allParameters;

    double centerMassX;
    double centerMassY;
    double centerMassdX;
    double centerMassdY;

    int nParticles;
    int nParticlesLong;
    int nParticlesTransv;
    int distrType;

    double channelRelRadius;
    double outEnergy;
    double voltage;
    double mass;
    double charge;
    double massNumber;
    double chargeNumber;

    double emittanceXunnorm;
    double emittanceYunnorm;

    double emittanceX;
    double emittanceY;

    double averagePhase;
    double phaseSpread;

    double Xmax;
    double dXmax;
    double GammaX;

    double Ymax;
    double dYmax;
    double GammaY;

    double impulseCurrent;

    double averageEnergy;
    double momentumSpread;

    double getMomentumSpread()
    {
        return momentumSpread;
    };

    double GetEmittanceX()
    {
        return emittanceX;
    };
    double GetEmittanceY()
    {
        return emittanceY;
    };
    double GetimpulseCurrent()
    {
        return impulseCurrent;
    };
    double GetchannelRelRadius()
    {
        return channelRelRadius;
    };
    double GetVoltage()
    {
        return voltage;
    };
    int GetnParticles()
    {
        return nParticles;
    };
    double GetCharge()
    {
        return charge;
    };
    double GetMass()
    {
        return mass;
    };
    void getTransvData(std::vector<double>& x, std::vector<double>& dx, std::vector<double>& y,
                       std::vector<double>& dy);
    void GetParametersAccelFlow(std::vector<std::string>& keysO, std::vector<double>& p);
    void SetParametersAccelFlow(const std::vector<double>& in);
    std::string getParticleType();
    arma::vec   GetTwissVector();
    double      GetTotalEnergy() const;
    double      getRestMassInGeV();
};
class SynchrotronFlow : public AccelFlowBase
{
  public:
    void GenerateParticlesAccel();
    void GenerateParticlesAccelTest(int nparticles, double xMax, double yMax, double dxMax,
                                    double dyMax);
    SynchrotronDynamics& GetDynamicsAccel()
    {
        return dynamicsAccel;
    };
    void GetMassCenterVector(std::vector<double>& result);
	void SetMassCenterVector(const std::vector<double>& result, const bool is_x);

    void SaveBetaToMadFile(const std::string& fileName);
    void SaveToMadFile(const std::string& fileName);
    SynchrotronFlow() : AccelFlowBase(0)
    {
    }
    SynchrotronDynamics& GetDynamicsAccelTest()
    {
        return dynamicsAccelTest;
    };

  private:
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    SynchrotronDynamics dynamicsAccel;
    SynchrotronDynamics dynamicsAccelTest;
};
class AccelFlow : public AccelFlowBase
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    LinacDynamicsLong   dynamicsLong;
    LinacDynamicsTransv dynamicsTransv;
    LinacDynamicsAccel  dynamicsAccel;
    LinacDynamicsAccel  dynamicsAccelTest;
    std::vector<double> oneParticle;
    std::vector<double> oneParticleData;

  public:
    LinacDynamicsTransv& GetdynamicsTransv()
    {
        return dynamicsTransv;
    };
    LinacDynamicsAccel& GetDynamicsAccel()
    {
        return dynamicsAccel;
    };
    LinacDynamicsAccel* GetDynamicsAccelPointer()
    {
        return &dynamicsAccel;
    };
    LinacDynamicsAccel& GetDynamicsAccelTest()
    {
        return dynamicsAccelTest;
    };

    AccelFlow(){};

    AccelFlow(int ProblemType) : AccelFlowBase(ProblemType){};

    void GenerateParticles(double lambda, double z0, double Rchann);
    void GenerateParticlesEnvelopes(double lambda, double z0, double Rchann);
    void GenerateParticlesAccel();

    double GetStartVelocity() const;
    double GetStartMomentum() const;
    double GetAlpha() const;
    double GetRestEnergy() const;
    double GetOutVelocity() const;
};

#endif