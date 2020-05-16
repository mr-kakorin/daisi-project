/*#ifndef LINACTOOLS_H
#define LINACTOOLS_H

#include "../common/Results.h"
#include "Dmath.h"
#include "LinacFlow.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <common_tools/constants.h>
#include <vector>
//#include "LinacModels.h"

namespace LinacTools
{
std::vector<double> dZtmp;

class RFQExternalGrid
{
  public:
    std::vector<double> Ez;
    std::vector<double> EzCol;
    std::vector<double> Ex;
    std::vector<double> Ey;
    std::vector<double> z;
    std::vector<double> rho;
    void CreateRFQGrid(const std::vector<double>& u1cells, const std::vector<double>& L,
                       const std::vector<double>& LL, const std::vector<double>& MinumalRadii,
                       int Ltotal, double lambda, double dz, double voltage);
};
void chargeCalculate(LinacDynamicsTransv& dyn, std::vector<double>& W1, std::vector<double>& rho,
                     double r);
void EzFieldSimulate(int centralCell, int ncells, RFQExternalGrid& grid, double Rb, double Rch);
void InCell(const std::vector<double>& z, LinacDynamicsTransv& dyn, std::vector<double>& W1,
            std::vector<double>& W2, std::vector<double>& rho);
void FieldInterpolate(LinacDynamicsTransv& dyn, RFQExternalGrid& grid, std::vector<double>& W1,
                      double time, double freq, double ph0, double I, int flag);
void calculateF(LinacDynamicsTransv& dyn, double alpha, std::vector<std::vector<double>>& f);
void calculateF(std::vector<double>& dyn, std::vector<double>& oneParticleData, double alpha,
                std::vector<double>& f);
void updateMomentumsAndPositions(std::vector<std::vector<double>>& dyn,
                                 std::vector<double>& oneParticleData, double dt, double alpha,
                                 int step);
void updateMomentumsAndPositionsCorrect(std::vector<std::vector<double>>& dyn,
                                        std::vector<double>& oneParticleData, double dt,
                                        double alpha, int step);

void InCell(const std::vector<double>& z, std::vector<double>& dyn,
            std::vector<double>& oneParticleData, RFQExternalGrid& grid, double time, double freq,
            double ph0, double I);
void InCell(const std::vector<double>& z, LinacDynamicsLong& dyn, std::vector<double>& W1,
            std::vector<double>& W2, std::vector<double>& rho);
void FieldInterpolate(LinacDynamicsLong& dyn, RFQExternalGrid& grid, std::vector<double>& W1,
                      double time, double freq, double ph0, double I);

void updateMomentums(LinacDynamicsLong& dyn, double dt, double alpha);
void updatePositions(LinacDynamicsLong& dyn, double dt, double alpha);
void updateMomentumsAndPositions(std::vector<LinacDynamicsTransv>& dyn, double dt, double alpha,
                                 int step, std::vector<std::vector<std::vector<double>>>& F);

void changeRFQcontrols(std::vector<std::vector<std::vector<double>>>& RFQControlsInput);
void GenerateRFQControlsbyParameterss(const std::vector<double>&       RFQParameters,
                                      std::vector<double>&             Modulations,
                                      std::vector<double>&             SyncPhases,
                                      std::vector<double>&             MatcherRadii,
                                      std::vector<double>&             OutputRadii,
                                      std::vector<std::vector<double>> RFQApproxParameters);

class LinacDynSimulator
{
    std::vector<double> W1;
    std::vector<double> W2;

    std::vector<double>                           parameters;
    std::vector<std::vector<std::vector<double>>> F;
    RFQExternalGrid                               externalGrid;

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& parameters;
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& parameters;
        parameters.resize(7);
    }

    void SetFieldSolverParameters();

  public:
    void SetParameters(std::vector<double> in)
    {
        parameters = in;
        SetFieldSolverParameters();
    };

    std::vector<double> GetParameters()
    {
        return parameters;
    };
    void CalculateRFQAcceptance(
        int flagType, LinacFlow& flow, const std::vector<double>& RFQParameters,
        const std::vector<double>& L, const std::vector<double>& Modulations,
        const std::vector<double>& MinumalRadii, const std::vector<double>& MatcherRadii,
        const std::vector<double>& AccEff, std::vector<DynamicsData*>& outputData);

    void SimulateRFQFlowDynamics(int flagType, int nCells, LinacFlow& flow,
                                 const std::vector<double>&  RFQParameters,
                                 const std::vector<double>&  L,
                                 const std::vector<double>&  MinumalRadii,
                                 const std::vector<double>&  AccEff,
                                 std::vector<DynamicsData*>& outputData, double startPhase);

    void SimulateRFQFlowDynamics(int startCell, int flagType, int nCells, LinacFlow& flow,
                                 const std::vector<double>&  RFQParameters,
                                 const std::vector<double>&  L,
                                 const std::vector<double>&  MinumalRadii,
                                 const std::vector<double>&  AccEff,
                                 std::vector<DynamicsData*>& outputData, double startPhase);

    void SimulateRFQLongDynamics(int startCell, int flagType, int nCells, LinacFlow& flow,
                                 const std::vector<double>&  RFQParameters,
                                 const std::vector<double>&  L,
                                 const std::vector<double>&  MinumalRadii,
                                 const std::vector<double>&  AccEff,
                                 std::vector<DynamicsData*>& outputData, double& Transmission,
                                 double startPhase);

    void SimulateOneParticleRFQ(int timeSignum, int startCell, int flagType, int nCells,
                                LinacFlow& flow, const std::vector<double>& RFQParameters,
                                const std::vector<double>& u1cells, const std::vector<double>& L,
                                const std::vector<double>& LL,
                                const std::vector<double>& MinumalRadii, int& result,
                                double startPhase, std::vector<DynamicsData*>& outputData);

    LinacDynSimulator();
    void Init();
};

class RFQFitnessDynamics
{
  public:
    LinacDynSimulator      dynamicsSimulator;
    std::vector<LinacFlow> flows;
    std::vector<double>    RFQParameters;
    std::vector<double>    OutRadii;

    std::vector<double> CellsLengths;
    std::vector<double> MinumalRadii;
    std::vector<double> MatcherRadii;
    std::vector<double> AccEff;

    std::vector<double> Modulations;
    std::vector<double> SyncPhases;

    std::vector<double>              AllRadii;
    std::vector<double>              AvEnergies;
    std::vector<DynamicsData*>       outputData;
    std::vector<std::vector<double>> RFQApproxParameters;

    RFQFitnessDynamics(){

    };
    RFQFitnessDynamics(LinacDynSimulator& dynamicsSimulatorIn, std::vector<LinacFlow>& flowsIn,
                       std::vector<double>&             RFQParametersIn,
                       std::vector<std::vector<double>> RFQApproxParametersIn,
                       std::vector<double>&             MatcherRadiiIn)
    {
        dynamicsSimulator   = dynamicsSimulatorIn;
        flows               = flowsIn;
        RFQParameters       = RFQParametersIn;
        RFQApproxParameters = RFQApproxParametersIn;
        MatcherRadii        = MatcherRadiiIn;
    };

    double FitnessFunction(const std::vector<double>& population);
    void GeneratePopulation(std::vector<double>& population);
    void CheckPopulation(std::vector<double>& population);
};

class RFQFitnessDynamicsEm
{
  public:
    LinacDynSimulator      dynamicsSimulator;
    std::vector<LinacFlow> flows;
    std::vector<double>    RFQParameters;
    std::vector<double>    OutRadii;

    std::vector<double> CellsLengths;
    std::vector<double> MinumalRadii;
    std::vector<double> MatcherRadii;
    std::vector<double> AccEff;

    std::vector<double> Modulations;
    std::vector<double> SyncPhases;

    std::vector<double>              AllRadii;
    std::vector<double>              AvEnergies;
    std::vector<DynamicsData*>       outputData;
    std::vector<std::vector<double>> RFQApproxParameters;

    RFQFitnessDynamicsEm(){

    };
    RFQFitnessDynamicsEm(LinacDynSimulator& dynamicsSimulatorIn, std::vector<LinacFlow>& flowsIn,
                         std::vector<double>&             RFQParametersIn,
                         std::vector<std::vector<double>> RFQApproxParametersIn,
                         std::vector<double>&             MatcherRadiiIn)
    {
        dynamicsSimulator   = dynamicsSimulatorIn;
        flows               = flowsIn;
        RFQParameters       = RFQParametersIn;
        RFQApproxParameters = RFQApproxParametersIn;
        MatcherRadii        = MatcherRadiiIn;
    };

    double FitnessFunction(const std::vector<double>& population);
    void GeneratePopulation(std::vector<double>& population);
    void CheckPopulation(std::vector<double>& population);
};

class RFQFitnessDynamicsAcc
{
  public:
    LinacDynSimulator      dynamicsSimulator;
    std::vector<LinacFlow> flows;
    std::vector<double>    RFQParameters;
    std::vector<double>    OutRadii;

    std::vector<double> CellsLengths;
    std::vector<double> MinumalRadii;
    std::vector<double> MatcherRadii;
    std::vector<double> AccEff;

    std::vector<double> Modulations;
    std::vector<double> SyncPhases;

    std::vector<double>              AllRadii;
    std::vector<double>              AvEnergies;
    std::vector<DynamicsData*>       outputData;
    std::vector<std::vector<double>> RFQApproxParameters;

    RFQFitnessDynamicsAcc(){

    };
    RFQFitnessDynamicsAcc(LinacDynSimulator& dynamicsSimulatorIn, std::vector<LinacFlow>& flowsIn,
                          std::vector<double>&             RFQParametersIn,
                          std::vector<std::vector<double>> RFQApproxParametersIn,
                          std::vector<double>&             MatcherRadiiIn)
    {
        dynamicsSimulator   = dynamicsSimulatorIn;
        flows               = flowsIn;
        RFQParameters       = RFQParametersIn;
        RFQApproxParameters = RFQApproxParametersIn;
        MatcherRadii        = MatcherRadiiIn;
    };

    double FitnessFunction(const std::vector<double>& population);
    void GeneratePopulation(std::vector<double>& population);
    void CheckPopulation(std::vector<double>& population);
};

class RFQFitnessDynamicsMatcher
{
  public:
    LinacDynSimulator      dynamicsSimulator;
    std::vector<LinacFlow> flows;
    std::vector<double>    RFQParameters;
    std::vector<double>    OutRadii;

    std::vector<double> CellsLengths;
    std::vector<double> MinumalRadii;
    std::vector<double> MatcherRadii;
    std::vector<double> AccEff;

    std::vector<double> Modulations;
    std::vector<double> SyncPhases;

    std::vector<double>              AllRadii;
    std::vector<double>              AvEnergies;
    std::vector<DynamicsData*>       outputData;
    std::vector<std::vector<double>> RFQApproxParameters;

    RFQFitnessDynamicsMatcher(){

    };
    RFQFitnessDynamicsMatcher(LinacDynSimulator&      dynamicsSimulatorIn,
                              std::vector<LinacFlow>& flowsIn, std::vector<double>& RFQParametersIn,
                              std::vector<std::vector<double>> RFQApproxParametersIn,
                              std::vector<double>&             MatcherRadiiIn)
    {
        dynamicsSimulator   = dynamicsSimulatorIn;
        flows               = flowsIn;
        RFQParameters       = RFQParametersIn;
        RFQApproxParameters = RFQApproxParametersIn;
        MatcherRadii        = MatcherRadiiIn;
    };

    double FitnessFunction(const std::vector<double>& population);
    void GeneratePopulation(std::vector<double>& population);
    void CheckPopulation(std::vector<double>& population);
};

double InitialFitnessFunction(const std::vector<double>& population, RFQFitnessDynamics* dynObj);
};
#endif*/