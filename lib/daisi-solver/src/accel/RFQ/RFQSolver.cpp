#include <numeric>

#include "../base/AccelFlow.h"
#include "../base/AccelToolsGeneral.h"
#include "RFQDevice.h"
#include "RFQSolver.h"
#include "RFQStrings.h"
#include "RFQTools.h"
#include "Results.h"
#include "Tools.h"

#include <common_tools/constants.h>
#include <notk/controller.hpp>

std::vector<std::string>              empty;
std::vector<std::vector<std::string>> empty1;

double calc_fit_em(const std::vector<std::vector<float>>& S)
{
    double Sav    = 0;
    double result = 0;

    for (int np = 0; np < S.size(); np++)
    {
        Sav = Sav + S[np].back();
    }
    Sav = Sav / S.size();

    for (int np = 0; np < S.size(); np++)
    {
        result = result + std::pow(Sav - S[np].back(), 2.0);
    }
    return result;
}

RFQSolver::RFQSolver(const std::string& dataFolderIn) : AccelSolverBase(dataFolderIn, 0)
{
    std::vector<std::string>              empty;
    std::vector<std::vector<std::string>> empty1;

    AddSolver(GenerateRFQForFlowName, empty1, empty, empty1);

    AddSolver(SimRFQForFlowName, empty1, empty, empty1);

    AddSolver(ExportRFQForFlowName,
              std::vector<std::vector<std::string>>{
                  {}, {"Electrodes"}, {}, {}, {"File with electodes"}, {}},
              empty, empty1);

    AddSolver(
        MatcherOptName,
        std::vector<std::vector<std::string>>{{""}, {}, {"*.json"}, {"*.json"}, {}, {"Config"}},
        empty, empty1);

    std::vector<std::string> tt1 = {"Required energy", "Required length", "Search delta"};

    std::vector<std::vector<std::string>> tt7{{"Transmission", "Effective emittance"},
                                              {"Use matcher", "Don't use matcher"}};

    AddSolver(
        OptName,
        std::vector<std::vector<std::string>>{{""}, {}, {"*.json"}, {"*.json"}, {}, {"Config"}},
        tt1, tt7);

    BaseInit(dataFolderIn);

    /*inputFileNames = RFQInputFiles;
    outpuFileNames = RFQOutputFiles;
    inputFileExt = RFQInputFilesExt;
    inputFileProps = RFQInputFilesProps;
    solversNames = RFQSolversNames;
    BaseInit(dataFolderIn, RFQSolverParameteresNames);*/
};
void RFQSolver::GenerateRFQForFlow(std::shared_ptr<RFQDevice> device, double& progress,
                                   bool&                                              flagAbort,
                                   std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                   std::string&                                       errorMsg)
{

    int flow = 0;
    if (device->GetLinacFlows().size() == 0)
    {
        errorMsg = "There are no flows";
        return;
    };
    if (flow >= device->GetLinacFlows().size())
    {
        errorMsg = "Incorrect flow number";
        return;
    };
    int succes;

    std::vector<double> OutRadii;

    GenerateRFQForFlowPr(
        succes, device->GetMainAccelParameters(), device->GetLinacFlows()[flow],
        device->GetSequencesOfParameters()[0], device->GetSequencesOfParameters()[1],
        device->GetSequencesOfParameters()[2], device->GetSequencesOfParametersCalc()[0],
        device->GetSequencesOfParametersCalc()[1], device->GetSequencesOfParametersCalc()[2],
        OutRadii, device->GetSequencesOfParametersCalc()[3],
        device->GetSequencesOfParametersCalc()[4]);
};
void RFQSolver::SimulateDynamics(std::shared_ptr<RFQDevice> device, double& progress,
                                 bool&                                              flagAbort,
                                 std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                 std::string& errorMsg, std::vector<std::string>& status)
{

    int flow = 0;
    if (device->GetLinacFlows().size() == 0)
    {
        errorMsg = "There are no flows";
        return;
    };
    if (flow >= device->GetLinacFlows().size())
    {
        errorMsg = "Incorrect flow number";
        return;
    };
    outputData.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(RFQResults, SimRFQForFlowName, 6, {1, 1, 1, 1, 1, 1})));

    int nCellMatcher = device->GetSequencesOfParameters()[2].size();
    int nCellRegular = device->GetSequencesOfParameters()[0].size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double lambda  = commtools::LIGHT_VELOCITY() / device->GetMainAccelParameters()[0];
    double steps   = 150.0;
    double dt      = lambda / 150.0;
    double dzSteps = 0.5 * steps;
    std::shared_ptr<ExternalAccelGrid> grid =
        std::shared_ptr<ExternalAccelGrid>(new ExternalAccelGrid());
    grid->CreateRFQGrid(device->GetSequencesOfParametersCalc()[3],
                        device->GetSequencesOfParametersCalc()[0],
                        device->GetSequencesOfParametersCalc()[1], nCells, lambda, dzSteps,
                        device->GetLinacFlows()[flow]->GetVoltage());
    int nsaveTr = 100;
    SimulateFlowDynamicsEnvelopes(
        grid, steps, 0, 1, device->GetLinacFlows()[flow], device->GetMainAccelParameters()[0],
        device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
        outputData.back(), nsaveTr, device->GetSequencesOfParameters()[1][0], progress);

    float acc =
        std::min(outputData.back()->dataAdd[2].back(), outputData.back()->dataAdd[3].back());
    float L = outputData.back()->dataAdd[5].back();

    status.push_back("Acceleration = " + std::to_string(acc));
    status.push_back("length = " + std::to_string(L));
};

void RFQSolver::Opt(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                    std::string& errorMsg, std::vector<std::string>& status)
{
 /*   std::vector<std::shared_ptr<SimulationDataAccel>> outputDatatmp_s;
    outputDatatmp_s.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(RFQResults, SimRFQForFlowName, 6, {1, 1, 1, 1, 1, 1})));

    std::vector<std::vector<std::string>> keysO;
    std::vector<std::vector<double>>      p;

    device->GetSectionParameters(keysO, p);
    int flow = 0;

    std::vector<double> x_left;
    std::vector<double> x_right;
    std::vector<double> x_0;

    double delta = solverParameters[OptName]->find("Search delta");

    //	bool use_matcher =
    int use_matcher = solverParametersFlags[OptName][1];
    int fit_type    = solverParametersFlags[OptName][0];

    int start_section = 0;
    if (use_matcher == 1)
    {
        start_section = 1;
    }

    for (int section = start_section; section < p.size(); section++)
    {
        for (int i = 0; i < p[section].size(); i++)
        {
            x_left.push_back(p[section][i] * (1.0 - delta));
            x_right.push_back(p[section][i] * (1.0 + delta));
            x_0.push_back(p[section][i]);
        }
    }

    int nCellMatcher = device->GetSequencesOfParameters()[2].size();
    int nCellRegular = device->GetSequencesOfParameters()[0].size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double lambda  = commtools::LIGHT_VELOCITY() / device->GetMainAccelParameters()[0];
    double steps   = 150.0;
    double dt      = lambda / 150.0;
    double dzSteps = 0.5 * steps;

    int nsaveTr = 100;

    auto insert_par = [&](const std::vector<double>& input) {
        size_t k = 0;
        for (int section = start_section; section < p.size(); section++)
        {
            std::vector<double> params;

            for (size_t i = 0; i < p[section].size(); i++)
            {
                params.push_back(input[k]);
                k++;
            }

            device->SetAccelSectionParameters(section, params);
        }
    };

    auto req_en  = solverParameters[OptName]->find("Required energy");
    auto req_len = solverParameters[OptName]->find("Required length");

    std::shared_ptr<ExternalAccelGrid> grid =
        std::shared_ptr<ExternalAccelGrid>(new ExternalAccelGrid());
    grid->CreateRFQGrid(device->GetSequencesOfParametersCalc()[3],
                        device->GetSequencesOfParametersCalc()[0],
                        device->GetSequencesOfParametersCalc()[1], nCells, lambda, dzSteps,
                        device->GetLinacFlows()[flow]->GetVoltage());

    SimulateFlowDynamicsEnvelopes(
        grid, steps, 0, 1, device->GetLinacFlows()[flow], device->GetMainAccelParameters()[0],
        device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
        outputDatatmp_s.back(), nsaveTr, device->GetSequencesOfParameters()[1][0], progress);

    auto req_acc = std::min(*std::min_element(outputDatatmp_s.back()->dataAdd[2].begin(),
                                              outputDatatmp_s.back()->dataAdd[2].end()),
                            *std::min_element(outputDatatmp_s.back()->dataAdd[3].begin(),
                                              outputDatatmp_s.back()->dataAdd[3].end()));

    auto fitness = [&](const std::vector<double>& input) -> double {
        std::vector<std::shared_ptr<SimulationDataAccel>> outputDatatmp;
        outputDatatmp.push_back(std::shared_ptr<SimulationDataAccel>(
            new SimulationDataAccel(RFQResults, SimRFQForFlowName, 6, {1, 1, 1, 1, 1, 1})));

        insert_par(input);
        device->CreateSequences();
        if (!device->checkSequence())
        {
            return 1;
        }

        int succes;

        std::vector<double> OutRadii;

        GenerateRFQForFlowPr(
            succes, device->GetMainAccelParameters(), device->GetLinacFlows()[flow],
            device->GetSequencesOfParameters()[0], device->GetSequencesOfParameters()[1],
            device->GetSequencesOfParameters()[2], device->GetSequencesOfParametersCalc()[0],
            device->GetSequencesOfParametersCalc()[1], device->GetSequencesOfParametersCalc()[2],
            OutRadii, device->GetSequencesOfParametersCalc()[3],
            device->GetSequencesOfParametersCalc()[4]);

        auto en = device->GetSequencesOfParametersCalc()[4].back();
        auto L  = std::accumulate(device->GetSequencesOfParameters()[0].begin(),
                                 device->GetSequencesOfParameters()[0].begin(), 0);

        if (en < req_en || L > req_len)
        {
            return 1;
        }

        int nCellMatcher = device->GetSequencesOfParameters()[2].size();
        int nCellRegular = device->GetSequencesOfParameters()[0].size() - nCellMatcher;
        int nCells       = nCellRegular + nCellMatcher;

        std::shared_ptr<ExternalAccelGrid> grid =
            std::shared_ptr<ExternalAccelGrid>(new ExternalAccelGrid());
        grid->CreateRFQGrid(device->GetSequencesOfParametersCalc()[3],
                            device->GetSequencesOfParametersCalc()[0],
                            device->GetSequencesOfParametersCalc()[1], nCells, lambda, dzSteps,
                            device->GetLinacFlows()[flow]->GetVoltage());

        SimulateFlowDynamicsEnvelopes(
            grid, steps, 0, 1, device->GetLinacFlows()[flow], device->GetMainAccelParameters()[0],
            device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
            outputDatatmp.back(), nsaveTr, device->GetSequencesOfParameters()[1][0], progress);

        float acc = std::min(*std::min_element(outputDatatmp.back()->dataAdd[2].begin(),
                                               outputDatatmp.back()->dataAdd[2].end()),
                             *std::min_element(outputDatatmp.back()->dataAdd[3].begin(),
                                               outputDatatmp.back()->dataAdd[3].end()));
        // status.push_back("Acceleration = " + std::to_string(acc));
        if (fit_type == 0)
        {
            return 1.0 - acc;
        }
        if (fit_type == 1)
        {
            if (acc < req_acc)
            {
                return 1.0;
            }
            return calc_fit_em(outputDatatmp.back()->data[3]) +
                   calc_fit_em(outputDatatmp.back()->data[6]);
        }
    };

    notk::NOTKController<double, double> NOTKController;

    NOTKController.add_problem_config(inputFileNames[OptName][0]);

    NOTKController.set_borders(x_left, x_right);
    NOTKController.set_x_0(x_0);

    NOTKController.set_fitness(fitness);

    if (!NOTKController.check_configuration())
    {
        errorMsg = "An error occured during notk library initialization. Please "
                   "see log files for additional info.";
        return;
    }

    auto result = NOTKController.process(flagAbort, status);
    insert_par(result->get_last_argument());*/
}

void RFQSolver::MatcherOpt(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                           std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                           std::string& errorMsg, std::vector<std::string>& status)
{
   /* std::vector<std::shared_ptr<SimulationDataAccel>> outputDatatmp_s;
    outputDatatmp_s.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(RFQResults, SimRFQForFlowName, 6, {1, 1, 1, 1, 1, 1})));

    int flow = 0;
    if (device->GetLinacFlows().size() == 0)
    {
        errorMsg = "There are no flows";
        return;
    };
    if (flow >= device->GetLinacFlows().size())
    {
        errorMsg = "Incorrect flow number";
        return;
    };

    double lambda  = commtools::LIGHT_VELOCITY() / device->GetMainAccelParameters()[0];
    double steps   = 150.0;
    double dt      = lambda / 150.0;
    double dzSteps = 0.5 * steps;

    int nsaveTr = 100;

    int nCellMatcher = device->GetSequencesOfParameters()[2].size();
    int nCellRegular = device->GetSequencesOfParameters()[0].size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    std::shared_ptr<ExternalAccelGrid> grid =
        std::shared_ptr<ExternalAccelGrid>(new ExternalAccelGrid());
    grid->CreateRFQGrid(device->GetSequencesOfParametersCalc()[3],
                        device->GetSequencesOfParametersCalc()[0],
                        device->GetSequencesOfParametersCalc()[1], nCells, lambda, dzSteps,
                        device->GetLinacFlows()[flow]->GetVoltage());

    SimulateFlowDynamicsEnvelopes(
        grid, steps, 0, 1, device->GetLinacFlows()[flow], device->GetMainAccelParameters()[0],
        device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
        outputDatatmp_s.back(), nsaveTr, device->GetSequencesOfParameters()[1][0], progress);

    auto req_acc = std::min(*std::min_element(outputDatatmp_s.back()->dataAdd[2].begin(),
                                              outputDatatmp_s.back()->dataAdd[2].end()),
                            *std::min_element(outputDatatmp_s.back()->dataAdd[3].begin(),
                                              outputDatatmp_s.back()->dataAdd[3].end()));

    auto fitness = [&](const std::vector<double>& radii) -> double {
        std::vector<std::shared_ptr<SimulationDataAccel>> outputDatatmp;
        outputDatatmp.push_back(std::shared_ptr<SimulationDataAccel>(
            new SimulationDataAccel(RFQResults, SimRFQForFlowName, 6, {1, 1, 1, 1, 1, 1})));

        for (size_t i = 1; i < radii.size(); i++)
        {
            if (0.9999 * radii[i] > radii[i - 1])
            {
                return 1;
            }
        }

        int succes;

        std::vector<double> OutRadii;

        GenerateRFQForFlowPr(
            succes, device->GetMainAccelParameters(), device->GetLinacFlows()[flow],
            device->GetSequencesOfParameters()[0], device->GetSequencesOfParameters()[1], radii,
            device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
            device->GetSequencesOfParametersCalc()[2], OutRadii,
            device->GetSequencesOfParametersCalc()[3], device->GetSequencesOfParametersCalc()[4]);

        int nCellMatcher = device->GetSequencesOfParameters()[2].size();
        int nCellRegular = device->GetSequencesOfParameters()[0].size() - nCellMatcher;
        int nCells       = nCellRegular + nCellMatcher;

        std::shared_ptr<ExternalAccelGrid> grid =
            std::shared_ptr<ExternalAccelGrid>(new ExternalAccelGrid());
        grid->CreateRFQGrid(device->GetSequencesOfParametersCalc()[3],
                            device->GetSequencesOfParametersCalc()[0],
                            device->GetSequencesOfParametersCalc()[1], nCells, lambda, dzSteps,
                            device->GetLinacFlows()[flow]->GetVoltage());

        SimulateFlowDynamicsEnvelopes(
            grid, steps, 0, 1, device->GetLinacFlows()[flow], device->GetMainAccelParameters()[0],
            device->GetSequencesOfParametersCalc()[0], device->GetSequencesOfParametersCalc()[1],
            outputDatatmp.back(), nsaveTr, device->GetSequencesOfParameters()[1][0], progress);

        float acc = std::min(*std::min_element(outputDatatmp.back()->dataAdd[2].begin(),
                                               outputDatatmp.back()->dataAdd[2].end()),
                             *std::min_element(outputDatatmp.back()->dataAdd[3].begin(),
                                               outputDatatmp.back()->dataAdd[3].end()));

        if (acc < req_acc)
        {
            return 1.0;
        }

        return calc_fit_em(outputDatatmp.back()->data[3]) +
               calc_fit_em(outputDatatmp.back()->data[6]);
    };

    notk::NOTKController<double, double> NOTKController;

    NOTKController.add_problem_config(inputFileNames[MatcherOptName][0]);

    std::vector<double> current = device->GetSequencesOfParameters()[2];
    std::vector<double> x_left(current.size());
    std::fill(x_left.begin(), x_left.end(), current.back());

    std::vector<double> x_right(current.size());
    std::fill(x_right.begin(), x_right.end(), current[0] * 1.5);

    NOTKController.set_borders(x_left, x_right);
    NOTKController.set_x_0(current);

    NOTKController.set_fitness(fitness);

    if (!NOTKController.check_configuration())
    {
        errorMsg = "An error occured during notk library initialization. Please "
                   "see log files for additional info.";
        return;
    }

    auto result                           = NOTKController.process(flagAbort, status);
    device->GetSequencesOfParameters()[2] = result->get_last_argument();*/
}

void RFQSolver::ExportToLidos(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                              std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                              std::string&                                       errorMsg)
{
    std::vector<double>& MatcherRadii = device->GetSequencesOfParameters()[2];
    std::vector<double>& SyncPhases   = device->GetSequencesOfParameters()[1];
    std::vector<double>& Modulations  = device->GetSequencesOfParameters()[0];
    double               lambda = commtools::LIGHT_VELOCITY() / device->GetMainAccelParameters()[0];

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCell        = nCellRegular + nCellMatcher;
    // RFQParameters[6] = std::abs(f2(RFQParameters[3], kOut, RFQParameters[1], 0));
    int succes;
    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations,
    // SyncPhases, CellsLengths, MinumalRadii, MinumalRadiiRegular, MatcherRadii, AccEff,
    // AvEnergies);

    std::vector<double> LTmp;
    std::vector<double> MinumalRadiiTmp;
    std::vector<double> MinumalRadiiRegularTmp;
    std::vector<double> MatcherRadiiTmp;
    std::vector<double> OutRadiiTmp;
    std::vector<double> AccEffTmp;

    std::vector<double> OutRadii;

    GenerateRFQForFlowPr(
        succes, device->GetMainAccelParameters(), device->GetLinacFlows()[0],
        device->GetSequencesOfParameters()[0], device->GetSequencesOfParameters()[1],
        device->GetSequencesOfParameters()[2], device->GetSequencesOfParametersCalc()[0],
        device->GetSequencesOfParametersCalc()[1], device->GetSequencesOfParametersCalc()[2],
        OutRadii, device->GetSequencesOfParametersCalc()[3],
        device->GetSequencesOfParametersCalc()[4]);

    std::vector<double>& CellsLengths  = device->GetSequencesOfParametersCalc()[0];
    std::vector<double>& MinumalRadii  = device->GetSequencesOfParametersCalc()[1];
    std::vector<double>& AccEff        = device->GetSequencesOfParametersCalc()[3];
    std::vector<double>& AvEnergiesTmp = device->GetSequencesOfParametersCalc()[4];

    FILE* fp = fopen("d:/Dropbox/myWorks/rfq/parmteq.dat", "w");

    fprintf(fp, "%.0f %.0f %.6f %.6f %.6f\n", device->GetLinacFlows()[0]->chargeNumber,
            device->GetLinacFlows()[0]->massNumber, lambda,
            device->GetLinacFlows()[0]->impulseCurrent * 1000.0,
            device->GetLinacFlows()[0]->emittanceX / 10.0);
    fprintf(fp, "%d %d\n", nCellRegular + nCellMatcher - 1, nCellMatcher);

    double lCurrent = 0;
    double m;
    double LCell;
    double F = 0.5;
    double vel;
    double energy;
    for (int i = 0; i < nCellRegular + nCellMatcher; i++)
    {
        energy = AvEnergiesTmp[i] * 1e-6;

        vel = sqrt(2 * AvEnergiesTmp[i] * std::abs(commtools::ELECTRON_CHARGE()) /
                   device->GetLinacFlows()[0]->mass) /
              commtools::LIGHT_VELOCITY();

        if (i == 0)
            LCell = 0;
        else
            LCell = CellsLengths[i - 1] * 100;
        fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", i, lCurrent,
                MinumalRadii[i] * 100, Modulations[i], energy, vel,
                device->GetLinacFlows()[0]->voltage / 1000, AccEff[i],
                SyncPhases[i] * 180 / commtools::PI(), F, LCell);
        lCurrent = lCurrent + LCell;
    };
    fclose(fp);
};
