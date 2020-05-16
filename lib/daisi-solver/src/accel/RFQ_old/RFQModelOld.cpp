#include "LinacModels.h"
#include "OptimizationTools.h"

double f2(double m, double k, double r, double t)
{
    return t - commtools::PI() / 4 * (m * m - 1) /
                   (m * m * boost::math::cyl_bessel_i(0, 2 * k * r / (m + 1)) +
                    boost::math::cyl_bessel_i(0, 2 * m * k * r / (m + 1)));
}
void RFQDesignModel::ExportToLidos(std::string filename, int flow)
{
    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCell        = nCellRegular + nCellMatcher;
    // RFQParameters[6] = std::abs(f2(RFQParameters[3], kOut, RFQParameters[1], 0));
    int succes;
    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations, SyncPhases, CellsLengths,
    // MinumalRadii, MinumalRadiiRegular, MatcherRadii, AccEff, AvEnergies);

    std::vector<double> LTmp;
    std::vector<double> MinumalRadiiTmp;
    std::vector<double> MinumalRadiiRegularTmp;
    std::vector<double> MatcherRadiiTmp;
    std::vector<double> OutRadiiTmp;
    std::vector<double> AccEffTmp;
    std::vector<double> AvEnergiesTmp;

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations, SyncPhases, LTmp,
                                   MinumalRadiiTmp, MinumalRadiiRegularTmp, MatcherRadiiTmp, OutRadiiTmp, AccEffTmp,
                                   AvEnergiesTmp);

    FILE* fp;
    fopen_s(&fp, filename.c_str(), "w");

    fprintf_s(fp, "%.0f %.0f %.6f %.6f %.6f\n", LinacFlows[flow].chargeNumber, LinacFlows[flow].massNumber, lambda,
              LinacFlows[flow].impulseCurrent * 1000.0, LinacFlows[flow].emittanceX / 10.0);
    fprintf_s(fp, "%d %d\n", nCellRegular + nCellMatcher - 1, nCellMatcher);

    double lCurrent = 0;
    double m;
    double LCell;
    double F = 0.5;
    double vel;
    double energy;
    for (int i = 0; i < nCellRegular + nCellMatcher; i++)
    {
        energy = AvEnergiesTmp[i] * 1e-6;

        vel =
            sqrt(2 * AvEnergiesTmp[i] * std::abs(commtools::ELECTRON_CHARGE()) / LinacFlows[flow].mass) / commtools::LIGHT_VELOCITY();

        if (i == 0)
            LCell = 0;
        else
            LCell = CellsLengths[i - 1] * 100;
        fprintf_s(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", i, lCurrent,
                  MinumalRadii[i] * 100, Modulations[i], energy, vel, LinacFlows[flow].voltage / 1000, AccEff[i],
                  SyncPhases[i] * 180 / commtools::PI(), F, LCell);
        lCurrent = lCurrent + LCell;
    };
    fclose(fp);
};
void RFQDesignModel::SearchStartParameters(bool& flagAbort, double& progress, std::vector<double> params)
{
    simulationData.dataFlags.resize(1);
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + "0";
    simulationData.AddFlags(names);

    simulationData.reset();

    ParticlesSwarmOptimization<LinacTools::RFQFitnessDynamics> solver;
    int                                                        numThreads = params[0];
    LinacTools::RFQFitnessDynamics fitnessDyn(DynSimulator, LinacFlows, RFQParameters, RFQApproxParameters,
                                              MatcherRadii);
    solver.init(fitnessDyn, numThreads, params[1]);
    std::vector<double> parameters = solver.SearchMin(flagAbort, simulationData);

    progress = 1;

    int k = 0;
    for (int i = 0; i < RFQApproxParameters.size() - 1; i++)
    {
        for (int j = 0; j < RFQApproxParameters[i].size(); j++)
        {
            RFQApproxParameters[i][j] = parameters[k];
            k++;
        }
    };

    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases, MatcherRadii, OutRadii,
                                                 RFQApproxParameters);
    int succes;

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[0], Modulations, SyncPhases, CellsLengths,
                                   MinumalRadii, MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff, AvEnergies);

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double Transmission;
    DynSimulator.SimulateRFQLongDynamics(0, 0, nCells, LinacFlows[0], RFQParameters, CellsLengths, MinumalRadii, AccEff,
                                         outputDataLinac, Transmission, SyncPhases[0]);
};
void RFQDesignModel::RFQAccOpt(bool& flagAbort, double& progress, std::vector<double> params)
{
    simulationData.dataFlags.resize(1);
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + "0";
    simulationData.AddFlags(names);

    simulationData.reset();

    ParticlesSwarmOptimization<LinacTools::RFQFitnessDynamicsAcc> solver;
    int                                                           numThreads = params[0];

    LinacTools::RFQFitnessDynamicsAcc fitnessDyn(DynSimulator, LinacFlows, RFQParameters, RFQApproxParameters,
                                                 MatcherRadii);
    solver.init(fitnessDyn, numThreads, params[1]);
    std::vector<double> parameters = solver.SearchMin(flagAbort, simulationData);

    progress = 1;

    int k = 0;
    for (int i = 0; i < RFQApproxParameters.size(); i++)
    {
        for (int j = 0; j < RFQApproxParameters[i].size(); j++)
        {
            RFQApproxParameters[i][j] = parameters[k];
            k++;
        }
    };

    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases, MatcherRadii, OutRadii,
                                                 RFQApproxParameters);
    int succes;

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[0], Modulations, SyncPhases, CellsLengths,
                                   MinumalRadii, MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff, AvEnergies);

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double Transmission;
    DynSimulator.SimulateRFQLongDynamics(0, 0, nCells, LinacFlows[0], RFQParameters, CellsLengths, MinumalRadii, AccEff,
                                         outputDataLinac, Transmission, SyncPhases[0]);
};

void RFQDesignModel::RFQMatcherOpt(bool& flagAbort, double& progress, std::vector<double> params)
{
    simulationData.dataFlags.resize(1);
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + "0";
    simulationData.AddFlags(names);

    simulationData.reset();

    ParticlesSwarmOptimization<LinacTools::RFQFitnessDynamicsMatcher> solver;
    int                                                               numThreads = params[0];

    LinacTools::RFQFitnessDynamicsMatcher fitnessDyn(DynSimulator, LinacFlows, RFQParameters, RFQApproxParameters,
                                                     MatcherRadii);
    solver.init(fitnessDyn, numThreads, params[1]);
    std::vector<double> parameters = solver.SearchMin(flagAbort, simulationData);

    progress = 1;

    MatcherRadii = parameters;

    int succes;

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[0], Modulations, SyncPhases, CellsLengths,
                                   MinumalRadii, MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff, AvEnergies);

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double Transmission;
    DynSimulator.SimulateRFQLongDynamics(0, 0, nCells, LinacFlows[0], RFQParameters, CellsLengths, MinumalRadii, AccEff,
                                         outputDataLinac, Transmission, SyncPhases[0]);
};

void RFQDesignModel::RFQEmOpt(bool& flagAbort, double& progress, std::vector<double> params)
{
    simulationData.dataFlags.resize(1);
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + "0";
    simulationData.AddFlags(names);

    simulationData.reset();

    ParticlesSwarmOptimization<LinacTools::RFQFitnessDynamicsEm> solver;
    int                                                          numThreads = params[0];

    LinacTools::RFQFitnessDynamicsEm fitnessDyn(DynSimulator, LinacFlows, RFQParameters, RFQApproxParameters,
                                                MatcherRadii);
    solver.init(fitnessDyn, numThreads, params[1]);
    std::vector<double> parameters = solver.SearchMin(flagAbort, simulationData);

    progress = 1;

    int k = 0;
    for (int i = 0; i < RFQApproxParameters.size(); i++)
    {
        for (int j = 0; j < RFQApproxParameters[i].size(); j++)
        {
            RFQApproxParameters[i][j] = parameters[k];
            k++;
        }
    };

    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases, MatcherRadii, OutRadii,
                                                 RFQApproxParameters);
    int succes;

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[0], Modulations, SyncPhases, CellsLengths,
                                   MinumalRadii, MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff, AvEnergies);

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double Transmission;
    DynSimulator.SimulateRFQLongDynamics(0, 0, nCells, LinacFlows[0], RFQParameters, CellsLengths, MinumalRadii, AccEff,
                                         outputDataLinac, Transmission, SyncPhases[0]);
};

void RFQDesignModel::CalculateRFQAcc(int flow)
{
    int succes;
    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations, SyncPhases, CellsLengths,
    // MinumalRadii, MinumalRadiiRegular, MatcherRadii, AccEff, AvEnergies);  if (succes == 0) 	return;

    DynSimulator.Init();

    DynSimulator.CalculateRFQAcceptance(1, LinacFlows[flow], RFQParameters, CellsLengths, Modulations, MinumalRadii,
                                        MatcherRadii, AccEff, outputDataLinac);
};

void RFQDesignModel::SimulateLinacFlowDynamics(int flow)
{
    int succes;
    // LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations, SyncPhases, CellsLengths,
    // MinumalRadii, MinumalRadiiRegular, MatcherRadii, AccEff, AvEnergies);  if (succes == 0) 	return;
    DynSimulator.Init();

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCell        = nCellRegular + nCellMatcher;

    DynSimulator.SimulateRFQFlowDynamics(1, nCell, LinacFlows[flow], RFQParameters, CellsLengths, MinumalRadii, AccEff,
                                         outputDataLinac, SyncPhases[0]);
}

void RFQDesignModel::GenerateRFQForFlow(int& succes, int flow)
{
    /*std::vector<std::vector<std::vector<double>>> RFQControls0 = RFQControls;
    int nc = RFQParameters[2];

    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];*/

    LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], Modulations, SyncPhases, CellsLengths,
                                   MinumalRadii, MinumalRadiiRegular, MatcherRadii, OutRadii, AccEff, AvEnergies);

    /*while (1)
    {
            int failPeriod = LinacTools::GenerateRFQForFlow(succes, RFQParameters, LinacFlows[flow], RFQControls,
    CellsLengths, MinumalRadii, u1c, u2c, MatcherRadii, AllRadii, AvEnergies);
            


            double En = AvEnergies.back();

            if (succes == 1)
                    return;

            if (succes == 1 && En>LinacFlows[flow].outEnergy)
            {
                    return;
            }
            RFQControls = RFQControls0;
            RFQParameters[2] = nc;
            LinacTools::changeRFQcontrols(RFQControls);
            RFQParameters[2] = RFQControls[0][0].back()-1;
            /*if (succes == 0)
            {
                    for (int i = 0; i < RFQControls[0][1].size();i++)
                    {
                            if (RFQControls[0][0][i] > failPeriod)
                            {
                                    for (int j = i-1; j < i + 1;j++)
                                            RFQControls[0][1][j] = RFQControls[0][1][j] * 0.9;

                                    break;
                            };
                    }
            };
    };*/
};
RFQDesignModel::RFQDesignModel()
{
    simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    RFQParameters.resize(2);
    RFQApproxParameters.resize(6);
    for (int i = 0; i < RFQApproxParameters.size(); i++)
        RFQApproxParameters[i].resize(5);

    RFQApproxParameters[4].resize(2);
    RFQApproxParameters[5].resize(2);
};

std::vector<double>& RFQDesignModel::GetControl(int flag)
{
    if (flag == 0)
        return Modulations;
    if (flag == 1)
        return SyncPhases;
};

void RFQDesignModel::SetControlFunction(std::string file, int flag){
    /*FILE* fp;

    std::string filenameTmp;
    char ss[250];
    std::vector<int> boundariesTmp;
    double tmp1, tmp2;
    fopen_s(&fp, file.c_str(), "r");

    int i = 0;

    RFQControls[flag].resize(2);

    while (fgets(ss, 250, fp))
    {
            sscanf_s(ss, "%lf%lf", &tmp1, &tmp2);
            RFQControls[flag][0].push_back(tmp1);
            RFQControls[flag][1].push_back(tmp2);

    }
    fclose(fp);*/
};

void RFQDesignModel::SetLinacParameters(std::vector<double> in, std::vector<std::vector<double>> inParams,
                                        std::string folder)
{
    /*0 "Frequency (Hz)"
    1 "Channel Radius",
    2 "Number of buncher cells",
    3 "Buncher phase",
    4 "Buncher modulation",
    5 "Number of former cells",
    6 "Former phase",
    7 "Former modulation",
    8 "Number of accel cells",
    9 "Matcher max radius",
    10 "Number of matcher cells"*/

    for (int i           = 0; i < in.size(); i++)
        RFQParameters[i] = in[i];

    RFQApproxParameters = inParams;

    /*std::vector<double> parameters(9);

    for (int i = 0; i < 7; i++)
            parameters[i] = RFQParameters[i + 2];

    parameters[7] = 8.0;
    parameters[8] = 3.0;*/

    LinacTools::GenerateRFQControlsbyParameterss(RFQParameters, Modulations, SyncPhases, MatcherRadii, OutRadii,
                                                 RFQApproxParameters);

    /*std::vector<double> accEffStart;

    std::vector<double> phasesStart;


    FILE* fp;

    std::string f1 = folder + "\\u1.txt";
    std::string f2 = folder + "\\u2.txt";

    char ss[1000];

    double tmp;
    fopen_s(&fp, f1.c_str(), "r");

    while (fscanf(fp, "%lf", &tmp) > 0) {
            accEffStart.push_back(tmp);
    }

    fclose(fp);

    fopen_s(&fp, f2.c_str(), "r");

    while (fscanf(fp, "%lf", &tmp) > 0) {
            phasesStart.push_back(tmp);
    }

    fclose(fp);


    RFQControls.resize(2);
    RFQControls[0].resize(2);
    RFQControls[1].resize(2);

    RFQControls[0][1] = accEffStart;
    RFQControls[1][1] = phasesStart;

    RFQControls[0][0].resize(accEffStart.size());
    RFQControls[1][0].resize(phasesStart.size());

    RFQControls[0][0][0] = 0;
    RFQControls[1][0][0] = 0;

    int ncell = RFQParameters[2];

    int dn1 = ncell / (accEffStart.size() - 1);
    int dn2 = ncell / (phasesStart.size() - 1);


    for (int i = 0; i < accEffStart.size() - 2; i++)
            RFQControls[0][0][i + 1] = RFQControls[0][0][i] + dn1;

    RFQControls[0][0][accEffStart.size() - 1] = ncell;


    for (int i = 0; i < phasesStart.size() - 2; i++)
            RFQControls[1][0][i + 1] = RFQControls[1][0][i] + dn2;

    RFQControls[1][0][phasesStart.size() - 1] = ncell;*/
};

std::vector<double> RFQDesignModel::GetLinacParameters()
{
    return RFQParameters;
};

void RFQDesignModel::AddLinacFlow()
{
    LinacFlows.push_back(LinacFlow());
};
int RFQDesignModel::GetNumberOfLinacFlows()
{
    return LinacFlows.size();
}

void RFQDesignModel::GenerateMesh(std::string meshParam, double& progress)
{
    solver.meshGenerator.MeshGenerate(meshParam, progress, deviceStatus.GetDomainBoundary(), deviceStatus.mesh, 0);
    deviceStatus.mesh.Convert2GridData(deviceStatus.GetGridData());
};

void RFQDesignModel::AddFlow(int ParticleType, int DistributionStyle, double massNumber, double chargeNumber)
{
    int                      n     = GetNumberParticlesFlows();
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + std::to_string(n);
    simulationData.AddFlags(names);
    deviceStatus.AddFlow(ParticleType, DistributionStyle, massNumber, chargeNumber);

    if (0 == DistributionStyle)
        solver.emissionCurrentSolver.addFlow(n);
};

void RFQDesignModel::AddFlow(int ParticleType, int DistributionStyle)
{
    int                      n     = GetNumberParticlesFlows();
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + std::to_string(n);
    simulationData.AddFlags(names);
    deviceStatus.AddFlow(ParticleType, DistributionStyle);

    if (0 == DistributionStyle)
        solver.emissionCurrentSolver.addFlow(n);
};