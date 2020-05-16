// #include "LinacTools.h"
#include <random>
std::default_random_engine                    generator;
std::normal_distribution<double>              distribution(0, 1);
std::vector<std::vector<std::vector<double>>> F;
std::vector<std::vector<double>>              Fn;
std::vector<std::vector<double>>              FnOnePartilce;
std::vector<double>                           rightV;
std::vector<double>                           Vtmp;

double LinacTools::RFQFitnessDynamicsAcc::FitnessFunction(const std::vector<double>& population)
{
    std::vector<double> MinumalRadiiRegular;
    int                 k = 0;
    for (int i = 0; i < RFQApproxParameters.size(); i++)
    {
        for (int j = 0; j < RFQApproxParameters[i].size(); j++)
        {
            RFQApproxParameters[i][j] = population[k];
            k++;
        }
    };
    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases,
                                                 MatcherRadii, OutRadii, RFQApproxParameters);
    int succes;

    double result = 0;

    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, flows[0], Modulations, SyncPhases,
    // CellsLengths, MinumalRadii,
    //                                MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff,
    //                                AvEnergies);
    for (int i = 0; i < flows.size(); i++)
    {
        int nCellMatcher = MatcherRadii.size();
        int nCellRegular = Modulations.size() - nCellMatcher;
        int nCells       = nCellRegular + nCellMatcher;

        dynamicsSimulator.CalculateRFQAcceptance(0, flows[i], RFQParameters, CellsLengths,
                                                 Modulations, MinumalRadii, MatcherRadii, AccEff,
                                                 outputData);
        //	outputData.clear();
        //	dynamicsSimulator.SimulateRFQFlowDynamics(0, nCells, flows[i], RFQParameters,
        // CellsLengths,
        // MinumalRadii, AccEff, outputData, SyncPhases[0]);

        result = result + (5 - flows[i].emittanceX);
    }
    return result;
};
void LinacTools::RFQFitnessDynamicsAcc::CheckPopulation(std::vector<double>& population)
{
    population[5]  = int(population[5]);
    population[10] = int(population[10]);
    population[15] = int(population[15]);

    std::vector<int> nc  = {0, 5, 10, 15, 20};
    std::vector<int> deg = {3, 4, 8, 9, 13, 14, 18};

    for (int i = 0; i < nc.size(); i++)
    {
        population[nc[i]] = int(population[nc[i]]);
        if (population[nc[i]] < 5)
            population[nc[i]] = 5;
    }

    for (int i = 0; i < deg.size(); i++)
    {
        if (population[deg[i]] < 1)
            population[deg[i]] = 1;
    }

    // if (population[11] < -0.5)

    if (population[1] < -1.57)
        population[1] = -1.57;

    if (population[2] < 1.02)
        population[2] = 1.02;

    // if (population[12] < 1.48)
    population[11] = RFQApproxParameters[2][1];
    population[12] = RFQApproxParameters[2][2];
    population[15] = RFQApproxParameters[3][0];

    if (population[6] < population[1])
        population[6] = population[1];

    if (population[7] < population[2])
        population[7] = population[2];

    population[16] = population[11];
    population[17] = population[12];

    if (population[21] > 0.05)
        population[21] = 0.05;

    population[22] = 0;
    population[23] = 0;
};

void LinacTools::RFQFitnessDynamicsAcc::GeneratePopulation(std::vector<double>& population)
{
    population.clear();

    for (int i = 0; i < RFQApproxParameters.size(); i++)
    {
        for (int j = 0; j < RFQApproxParameters[i].size(); j++)
            population.push_back(RFQApproxParameters[i][j]);
    }

    /*std::vector<double> sigmas = { 0, 0, 0.05, 0, 0.5,
            0, 0, 0.05, 0, 0.5,
            0, 0, 0.05, 0, 0.5,
            0, 0, 0.05, 0, 0,
            1, 0.005, 0, 0 };*/

    std::vector<double> sigmas = {20,
                                  (population[1] + commtools::PI() / 2) / 3,
                                  (population[2] - 1) / 3,
                                  0.5,
                                  0.5,
                                  20,
                                  (population[6] + commtools::PI() / 2) / 3,
                                  (population[7] - 1) / 3,
                                  0.5,
                                  0.5,
                                  20,
                                  0.3,
                                  0.2,
                                  1,
                                  1.5,
                                  20,
                                  0,
                                  0,
                                  0,
                                  0};

    for (int i        = 0; i < population.size(); i++)
        population[i] = sigmas[i] * distribution(generator) + population[i];

    std::vector<int> nc  = {0, 5, 10, 15};
    std::vector<int> deg = {3, 4, 8, 9, 13, 14, 18, 19};

    for (int i = 0; i < nc.size(); i++)
    {
        population[nc[i]] = int(population[nc[i]]);
        if (population[nc[i]] < 5)
            population[nc[i]] = 5;
    }

    for (int i = 0; i < deg.size(); i++)
    {
        if (population[deg[i]] < 1)
            population[deg[i]] = 1;
    }
};
