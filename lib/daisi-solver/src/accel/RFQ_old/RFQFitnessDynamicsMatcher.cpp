#include "LinacTools.h"
#include <random>
std::default_random_engine                    generator;
std::normal_distribution<double>              distribution(0, 1);
std::vector<std::vector<std::vector<double>>> F;
std::vector<std::vector<double>>              Fn;
std::vector<std::vector<double>>              FnOnePartilce;
std::vector<double>                           rightV;
std::vector<double>                           Vtmp;

double LinacTools::RFQFitnessDynamicsMatcher::FitnessFunction(const std::vector<double>& population)
{
    int                 succes;
    std::vector<double> MinumalRadiiRegular;

    double result = 0;

    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases,
                                                 MatcherRadii, OutRadii, RFQApproxParameters);
    double lR           = MatcherRadii.back();
    MatcherRadii        = population;
    MatcherRadii.back() = lR;

    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, flows[0], Modulations, SyncPhases,
    // CellsLengths, MinumalRadii,
    //                                MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff,
    //                                AvEnergies);

    for (int i = 0; i < 1; i++)
    {
        int nCellMatcher = MatcherRadii.size();
        int nCellRegular = Modulations.size() - nCellMatcher;
        int nCells       = nCellRegular + nCellMatcher;

        outputData.clear();
        dynamicsSimulator.SimulateRFQFlowDynamics(0, nCells, flows[i], RFQParameters, CellsLengths,
                                                  MinumalRadii, AccEff, outputData, SyncPhases[0]);

        double Transmission = 1;
        for (int j = 0; j < outputData[0]->dataAdd[3].size(); j++)
        {
            if (outputData[0]->dataAdd[3][j] < Transmission)
                Transmission = outputData[0]->dataAdd[3][j];

            if (outputData[0]->dataAdd[2][j] < Transmission)
                Transmission = outputData[0]->dataAdd[2][j];
        };
        //	double Transmission = std::min(outputData[0]->dataAdd[3].back(),
        //outputData[0]->dataAdd[2].back());
        result = result + (1 - Transmission);
    };

    return result;
};
void LinacTools::RFQFitnessDynamicsMatcher::CheckPopulation(std::vector<double>& population)
{
    for (int i = 0; i < population.size() - 1; i++)
    {
        if (population[i + 1] > population[i])
        {
            population[i + 1] = population[i];
        }
    }
    /*population[5] = int(population[5]);
    population[10] = int(population[10]);
    population[15] = int(population[15]);

    std::vector<int> nc = { 0, 5, 10, 15, 20 };
    std::vector<int> deg = { 3, 4, 8, 9, 13, 14, 18 };

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

    //if (population[11] < -0.5)

    if (population[1] < -1.57)
            population[1] = -1.57;

    if (population[2] < 1.02)
            population[2] = 1.02;

    //if (population[12] < 1.48)
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
    population[23] = 0;*/
};

void LinacTools::RFQFitnessDynamicsMatcher::GeneratePopulation(std::vector<double>& population)
{

    population.clear();

    double mD = 2;

    double Rmax         = RFQApproxParameters[4][1];
    int    nCellMatcher = int(RFQApproxParameters[4][0]);

    double kMatcher = (Rmax - RFQParameters[1]) / (pow(nCellMatcher, mD));

    population.resize(nCellMatcher);

    for (int i = 0; i < nCellMatcher; i++)
    {
        population[i] = RFQParameters[1] + pow(nCellMatcher - i, mD) * kMatcher;
    }

    for (int i = 0; i < population.size() - 1; i++)
        population[i] =
            std::abs((population[i + 1] - population[i])) * distribution(generator) + population[i];
};
