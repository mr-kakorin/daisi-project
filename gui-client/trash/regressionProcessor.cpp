#include "regressionProcessor.h"
#include "unordered_map"
#include "vector"
#include <conio.h>
#include <exception>
#include <functional>
#include <iostream>
#include <random>
#include <string>

std::default_random_engine             generator;
std::normal_distribution<double>       distribution(0, 1);
std::uniform_real_distribution<double> distributionUn(0, 1);

template <class T> T vectormin(const std::vector<T>& v)
{
    float res = 20e20;
    for (int i = 0; i < v.size(); i++)
        if (v[i] < res)
            res = v[i];

    return res;
};

template <class T> T vectormax(const std::vector<T>& v)
{
    float res = -20e20;
    for (int i = 0; i < v.size(); i++)
        if (v[i] > res)
            res = v[i];

    return res;
};
void strsplit(char* string, char* split, std::vector<std::string>& result)
{
    result.clear();
    char* pch = strtok(string, split);

    while (pch != NULL)
    {
        result.push_back(pch);
        pch = strtok(NULL, split);
    }
};


void loadXls(std::string& errorMessage, std::unordered_map<std::string, std::vector<std::vector<float>>>& Oilwells)
{

    /// load data from csv file
    FILE*   fp;
    char    line[1250];
    errno_t err;

    err = fopen_s(&fp, "table.csv", "r");

    if (err != 0)
    {
        errorMessage = "file not found";
        return;
    };
    std::vector<std::string> tmpStr;
    std::vector<float>       tmp(3);
    fgets(line, 1250, fp);
    int ii = 0;
    while (fgets(line, 1250, fp))
    {
        ii++;
        strsplit(line, ",", tmpStr);
        if (!tmpStr[1].size())
            continue;

        auto search = Oilwells.find(tmpStr[1]);

        if (search == Oilwells.end())
        {
            search =
                Oilwells.insert(std::begin(Oilwells), std::make_pair(tmpStr[1], std::vector<std::vector<float>>(3)));
            // some preallocation for fast vectors filling
            for (int par = 0; par < 3; par++)
                search->second[par].reserve(10000);

            search->second[0].push_back(0);
        }

        bool flagAdd = true;
        // check for nonzero
        for (int par = 0; par < 3; par++)
        {
            tmp[par] = std::stof(tmpStr[par + 2]);
            if (0 == tmp[par])
            {
                flagAdd = false;
                break;
            }
        }
        if (flagAdd)
        {
            // we convert mounth volume to the mounth debit
            search->second[1].push_back(tmp[1] / tmp[0]);
            search->second[2].push_back(tmp[2] / tmp[0]);
            // we save time point where the debit is calculated
            search->second[0].push_back(search->second[0].back() + tmp[0]);
        };
    }
    fclose(fp);
};
void regressionFunctions(const std::vector<float>& xArray, std::vector<float>& yArray,
                         const std::vector<float>& parameters, int flag)
{
    /// regressionFunctions

    yArray.resize(xArray.size());
    switch (flag)
    {
    case 0:
        for (int i    = 0; i < xArray.size(); i++)
            yArray[i] = exp(-parameters[0] * xArray[i]);
        break;
    case 1:
        for (int i    = 0; i < xArray.size(); i++)
            yArray[i] = pow(double(1.0 + xArray[i]), double(-parameters[0]));
        break;
    case 2:
        for (int i = 0; i < xArray.size(); i++)
        {
            double t1 = xArray[i] * parameters[0] * parameters[1];
            double t2 = -1.0 / parameters[1];
            yArray[i] = pow(1.0 + t1, t2);
        }
        break;
    case 3:
        for (int i = 0; i < xArray.size(); i++)
        {
            if (xArray[i] <= parameters[0])
                yArray[i] = pow(double(1.0 + xArray[i]), -parameters[1]);
            else
                yArray[i] = pow(double(1.0 + (xArray[i] - parameters[0]) * parameters[1] * parameters[2]),
                                double(-1.0 / parameters[1]));
        }
        break;
    };
}
double regressionFitness(const std::unordered_map<std::string, std::vector<std::vector<float>>>& Oilwells, int flag,
                         int flagParameter, const std::vector<float>& parameters)
{
    // fitness function for regression
    double             result = 0;
    std::vector<float> yArray;
    for (auto it = Oilwells.begin(); it != Oilwells.end(); it++)
    {
        regressionFunctions(it->second[0], yArray, parameters, flag);
        for (int i = 0; i < it->second[1].size(); i++)
        {
            double tmp = log(it->second[flagParameter][i]) - log(it->second[flagParameter][0]) - log(yArray[i]);
            result     = result + tmp * tmp;
        }
    }
    return result;
};

void regressionEps(std::vector<float>& eps,
                   const std::unordered_map<std::string, std::vector<std::vector<float>>>& Oilwells, int flag,
                   int flagParameter, const std::vector<float>& parameters)
{
    // fitness function for regression
    double             result = 0;
    std::vector<float> yArray;
    for (auto it = Oilwells.begin(); it != Oilwells.end(); it++)
    {
        regressionFunctions(it->second[0], yArray, parameters, flag);
        for (int i = 0; i < it->second[1].size(); i++)
            eps.push_back(log(it->second[flagParameter][i]) - log(it->second[flagParameter][0]) - log(yArray[i]));
    }
};

template <class T>
void MMCMinimizationFunction(std::function<double(const std::vector<T>&)> fitnessFunction, std::vector<T> leftBorders,
                             std::vector<T> rightBorders, int maxIterations, std::vector<T>& result)
{
    // Monte-Carlo optimization routine

    double         bestValue = 1e20;
    std::vector<T> currentParameters(leftBorders.size());

    std::vector<T> sizes(leftBorders.size());
    for (int p   = 0; p < leftBorders.size(); p++)
        sizes[p] = rightBorders[p] - leftBorders[p];

    result.resize(leftBorders.size());
    for (int p    = 0; p < leftBorders.size(); p++)
        result[p] = leftBorders[p] + distributionUn(generator) * (rightBorders[p] - leftBorders[p]);

    for (int i = 0; i < maxIterations; i++)
    {

        // generate random values of parameters
        for (int p               = 0; p < leftBorders.size(); p++)
            currentParameters[p] = leftBorders[p] + distributionUn(generator) * (rightBorders[p] - leftBorders[p]);

        // fitness Function

        double currentFitness = fitnessFunction(currentParameters);
        if (currentFitness < bestValue)
        {
            result    = currentParameters;
            bestValue = currentFitness;
        };
        // contraction of intervals

        if ((i + 1) % 10 == 0)
        {
            for (int p = 0; p < leftBorders.size(); p++)
            {
                if (std::abs(result[p] - leftBorders[p]) < std::abs(result[p] - rightBorders[p]))
                    leftBorders[p] = leftBorders[p] + (rightBorders[p] - leftBorders[p]) * 0.2;
                else
                    rightBorders[p] = rightBorders[p] - (rightBorders[p] - leftBorders[p]) * 0.2;
            };
        };
        if (rightBorders[0] - leftBorders[0] < 1e-6 * sizes[0])
            break;
    };
};
void processOilwells(std::string& errorMessage, double& progress, std::vector<std::vector<float>>& results,
                     std::vector<std::vector<float>>& resultsLineX, std::vector<std::vector<float>>& resultsLineY,
                     std::vector<std::vector<float>>& eps)
{
    progress = 0.01;
    // we use std map with hashing for fast table loading
    std::unordered_map<std::string, std::vector<std::vector<float>>> Oilwells;
    loadXls(errorMessage, Oilwells);

    // collecting normalized statistical data for visualization
    results.resize(4);
    for (auto it = Oilwells.begin(); it != Oilwells.end(); it++)
    {
        for (int i = 0; i < it->second[1].size(); i++)
        {
            results[0].push_back(it->second[0][i]);
            ;
            results[1].push_back(it->second[1][i] / it->second[1][0]);
            results[2].push_back(it->second[2][i] / it->second[2][0]);
        }
    };

    // search parameters
    std::vector<std::vector<float>> result(4);

    MMCMinimizationFunction<float>(std::bind(regressionFitness, Oilwells, 0, 1, std::placeholders::_1), {0},
                                   {float(1e-5)}, 2e2, result[0]);
    MMCMinimizationFunction<float>(std::bind(regressionFitness, Oilwells, 2, 1, std::placeholders::_1), {0, 0},
                                   {float(1e-5), float(1e-5)}, 2e2, result[1]);

    MMCMinimizationFunction<float>(std::bind(regressionFitness, Oilwells, 0, 2, std::placeholders::_1), {0},
                                   {float(1e-5)}, 2e2, result[2]);
    MMCMinimizationFunction<float>(std::bind(regressionFitness, Oilwells, 1, 2, std::placeholders::_1), {0},
                                   {float(1e-5)}, 2e2, result[3]);

    // calculation of approximation functions
    float min      = vectormin(results[0]);
    float max      = vectormax(results[0]);
    int   NpointsV = 1000;
    float dx       = (max - min) / (NpointsV - 1);
    resultsLineX.resize(4);
    resultsLineY.resize(4);

    for (int plot = 0; plot < resultsLineX.size(); plot++)
        for (int i = 0; i < NpointsV; i++)
            resultsLineX[plot].push_back(min + i * dx);

    std::vector<int> plots = {0, 2, 0, 1};
    for (int plot = 0; plot < resultsLineX.size(); plot++)
        regressionFunctions(resultsLineX[plot], resultsLineY[plot], result[plot], plots[plot]);

    progress = 1.0;

    eps.resize(2);
    regressionEps(eps[0], Oilwells, 0, 1, result[0]);
    regressionEps(eps[1], Oilwells, 0, 2, result[2]);
};