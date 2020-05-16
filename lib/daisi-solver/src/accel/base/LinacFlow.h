#ifndef LinacFlow_H
#define LinacFlow_H

#include <vector>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

#include <common_tools/constants.h>
// class LinacDynamicsLong
// {
//   public:
//     int    NparticlesLong;
//     double longQ;

//     std::vector<double> p;
//     std::vector<double> z;
//     std::vector<int>    cellNumber;
//     std::vector<double> Ez;
//     void clearandResize(int nParticles);
//     std::vector<void*> GetData();
//     int checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
//                           const std::vector<double>& LL, const std::vector<double>& MinumalRadii,
//                           double lambda);
// };

// class LinacDynamicsTransv
// {
//   public:
//     int                 NparticlesTransv;
//     double              transvQ;
//     std::vector<double> p;
//     std::vector<double> z;
//     std::vector<double> S11x;
//     std::vector<double> S12x;
//     std::vector<double> S22x;
//     std::vector<double> S11y;
//     std::vector<double> S12y;
//     std::vector<double> S22y;
//     std::vector<int>    cellNumber;
//     std::vector<double> Ex;
//     std::vector<double> Ey;
//     std::vector<double> Ez;
//     std::vector<double> EzCol;

//     void clearandResize(int nParticles);
//     std::vector<void*> GetData();
//     int checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
//                           const std::vector<double>& LL, const std::vector<double>& MinumalRadii,
//                           double lambda, double relChann);
// };

class LinacFlow
{
    void translateParameters()
    {
        allParameters.resize(14);
        voltage      = allParameters[0];
        massNumber   = allParameters[1];
        chargeNumber = allParameters[2];
        mass         = massNumber * commtools::NEUTRON_MASS();
        charge       = -chargeNumber * commtools::ELECTRON_CHARGE();

        emittanceX     = allParameters[3];
        emittanceY     = allParameters[4];
        averageEnergy  = allParameters[5];
        momentumSpread = allParameters[6];
        impulseCurrent = allParameters[7];

        nParticlesLong   = int(allParameters[8]);
        nParticlesTransv = int(allParameters[9]);

        Xmax             = allParameters[10];
        dXmax            = allParameters[11];
        outEnergy        = allParameters[12];
        channelRelRadius = allParameters[13];
    };
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& allParameters;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& allParameters;
        translateParameters();
    }

  public:
    void GenerateParticlesRFQ(double lambda, double z0, double Rchann);
    void GenerateOneParticleRFQ(double lambda, double z0);

    double GetStartVelocity() const;
    double GetStartMomentum() const;
    double GetAlpha() const;
    double GetRestEnergy() const;
    double GetOutVelocity() const;

    std::vector<double> allParameters;

    int    nParticlesLong;
    int    nParticlesTransv;
    double channelRelRadius;
    double outEnergy;
    double voltage;
    double mass;
    double charge;
    double massNumber;
    double chargeNumber;

    double emittanceX;
    double Xmax;
    double dXmax;
    double GammaX;

    double emittanceY;
    double Ymax;
    double dYmax;
    double GammaY;

    double impulseCurrent;

    double averageEnergy;
    double momentumSpread;
    // LinacDynamicsLong   dynamicsLong;
    // LinacDynamicsTransv dynamicsTransv;
    std::vector<double> oneParticle;
    std::vector<double> oneParticleData;

    LinacFlow()
    {
        allParameters.resize(20);
    }

    std::vector<double> GetParametersRFQ()
    {
        return allParameters;
    };
    std::vector<double> GetParametersAPF()
    {
        return allParameters;
    };
    void SetParametersRFQ(std::vector<double> in)
    {
        allParameters = in;
        translateParameters();
    };
    void SetParametersAPF(std::vector<double> in)
    {
        allParameters = in;
        translateParameters();
    };
};

#endif