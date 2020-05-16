#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>
#include <daisi-solver/FlagStringsSolver.h>
#include <daisi-solver/SynchrotronStrings.h>


#include "VisualizationFunctionsAccel.h"
#include "VisualizationFunctionsNucl.h"
#include "colors.h"
#include "vtkComponent.h"

#include <memory>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include "VisualizationFunctionsAccel.h"
#include "VisualizationFunctionsNucl.h"
#include "colors.h"
#include "vtkComponent.h"

const std::vector<std::string> LinacResults = {
    "Energy(z)", "Phase(t)", "Rx(z)", "Ry(z)", "XdX Ellipses", "YdY Ellipses", "Acceleration and Transmission"};

const double PI                  = 3.141592653;
const double LIGHT_VELOCITY      = 299792458.0;
const double ELECTRON_CHARGE     = -1.60217656e-19;
const double ELECTRON_MASS       = 9.10938215e-31;
const double PROTON_MASS         = 1.67262178e-27;
const double VACUUM_PERMITTIVITY = 8.85418781e-12;
const double VACUUM_PERMEABILITY = 1.25663706e-6;
void         ShowLinacSequence(vtkComponent* input, AccelModelInterface* currentModel, int flag)
{
    input->CheckType(1);
    input->clear();
    float               ymax;
    float               ymin;
    std::vector<double> data;
    currentModel->GetSequencesOfParameters(flag, data);

    //	data = currentModel->GetControl(flag);

    std::vector<float> dataX(data.size());
    std::vector<float> dataY(data.size());

    for (int i = 0; i < data.size(); i++)
    {
        dataX[i] = i;
        dataY[i] = data[i];
    }

    std::string yName  = currentModel->GetnamesOfSequences()[flag];
    std::string xName  = "cell number";
    std::string contor = "Sequence";

    input->addVisualizationDataPlot(dataX, dataY, 2, xName, yName, contor);
};

void ShowAccelCalcParameters(vtkComponent* input, AccelModelInterface* currentModel, int flag)
{
    input->CheckType(1);
    input->clear();
    float               ymax;
    float               ymin;
    std::vector<double> data;
    currentModel->GetSequencesOfParametersCalc(flag, data);

    //	data = currentModel->GetControl(flag);

    std::vector<float> dataX(data.size());
    std::vector<float> dataY(data.size());

    for (int i = 0; i < data.size(); i++)
    {
        dataX[i] = i;
        dataY[i] = data[i];
    }

    std::string yName  = currentModel->GetcalcParametersNames()[flag];
    std::string xName  = "cell number";
    std::string contor = "Sequence";

    input->addVisualizationDataPlot(dataX, dataY, 2, xName, yName, contor);
};

void ShowSimulationNuclDataPlot(vtkComponent* input, AccelModelInterface* currentModel,
                                std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                                std::vector<int> props, std::vector<double> plotProperties)
{

    float ymax = -2e20;
    float ymin = 2e20;
    input->CheckType(1);
    input->clear();

    // int i = 0;
    ////input->addVisualizationDataFillRect({ 0.0,100.0 }, { 0.0,100.0 }, color, names[i], 0.5);

    int flagSeq = 0;

    if (simData->tag == SynchrotronSolversNameTwiss || simData->tag == SynchrotronSolversNameTwiss + " MADX")
        ShowSimulationNuclTwiss(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversBeam || simData->tag == SynchrotronSolversBeam + " MADX")
        ShowSimulationNuclBeam(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if ((simData->tag == SynchrotronSolversSVD) || (simData->tag == SynchrotronSolversMikado))
        ShowNuclCorrOrbit(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversNameCenterOfMass ||
        simData->tag == SynchrotronSolversNameCenterOfMass + " MADX")
        ShowSimulationNuclCM(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversOrbitCalc)
        ShowSimulationNuclOrbit(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversAcceptance)
        ShowSimulationNuclAcceptance(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversTolerancesNameReconstruction)
        SynchrotronRec(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (simData->tag == SynchrotronSolversCorrEffName)
        SynchrotronEst(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

	if (simData->tag == SynchrotronSolversMagnetShuffleName)
		SynchrotronShuffle(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);

    if (!flagSeq)
        ShowNuclElements(input, currentModel, simData, flag, props, plotProperties, flagSeq, ymax, ymin);
};

void ShowSimulationAccelDataPlot(vtkComponent* input, AccelModelInterface* currentModel,
                                 std::shared_ptr<SimulationDataAccel> data, const std::string& flag,
                                 std::vector<int> props, std::vector<double> plotProperties)
{
    input->CheckType(1);
    input->clear();
    std::string                     yName, xName;
    std::vector<std::vector<float>> tmp;
    float                           ymax;
    float                           ymin;
    float                           gamma;
    float                           en;
    std::vector<float>              Zav(data->data[1][0].size());
    std::vector<float>              beta           = data->dataAdd[4];
    int                             numberOfTraces = plotProperties[0];
    int                             Distrtype      = 0;
    if (flag == LinacResults[0])
    {
        yName = "Energy, Ev";
        xName = "Z, m";
        tmp   = data->data[1];

        en = LIGHT_VELOCITY * LIGHT_VELOCITY * data->mass;

        for (int i = 0; i < tmp.size(); i++)
        {
            std::vector<float> Xtmp;
            for (int j = 0; j < tmp[i].size(); j++)
            {
                float V   = tmp[i][j] * LIGHT_VELOCITY;
                tmp[i][j] = data->mass * V * V / (std::abs(ELECTRON_CHARGE) * 2);
            }
        }

        //	input->addVisualizationDataPlots(data->data[0], tmp, colors::blackcolor, 1, xName, yName,
        // numberOfTraces, Distrtype, props);
        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName,
                                              plotProperties[0], 0, props, ymax, ymin);
        return;
    }

    if (flag == LinacResults[1])
    {
        yName = "Phase, rad";
        xName = "T, ns";

        tmp = data->data[0];

        for (int k = 0; k < data->data[0][0].size(); k++)
        {
            Zav[k] = data->data[0][0][k];
            /*int s = 0;
            for (int i = 0; i < data->TimeArray.size(); i++)
            {
            if (k < data->data[0][i].size())
            {
            Zav[k] = Zav[k] + data->data[0][i][k];
            s++;
            }
            }
            Zav[k] = Zav[k] / s;*/

            for (int i = 0; i < data->TimeArray.size(); i++)
            {
                tmp[i][k] = tmp[i][k] - Zav[k];
                tmp[i][k] = 2 * PI * tmp[i][k] / (beta[k] * data->lambda);
                while (tmp[i][k] < -PI)
                    tmp[i][k] = tmp[i][k] + 2 * PI;
            }
        }

        input->addVisualizationDataPlots(data->TimeArray, tmp, colors::blackcolor, 1, xName, yName, plotProperties[0],
                                         0, props, -3.5, 3.5);
        return;
    };
    if (flag == LinacResults[2])
    {

        yName = "Rx, m";
        xName = "Z, m";
        tmp.resize(data->data[2].size());

        for (int i = 0; i < data->data[2].size(); i++)
        {
            tmp[i].resize(data->data[2][i].size());
            for (int j = 0; j < data->data[2][i].size(); j++)
                tmp[i][j] = sqrt(data->data[2][i][j]);
        }

        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                              Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlot(data->dataAdd[5], data->dataAdd[6], 2, xName, yName, std::string("Channel"));
        return;
    }
    if (flag == LinacResults[3])
    {
        yName = "Ry, m";
        xName = "Z, m";
        tmp.resize(data->data[5].size());
        for (int i = 0; i < data->data[5].size(); i++)
        {
            tmp[i].resize(data->data[5][i].size());
            for (int j = 0; j < data->data[5][i].size(); j++)
                tmp[i][j] = sqrt(data->data[5][i][j]);
        }

        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                              Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlot(data->dataAdd[5], data->dataAdd[6], 2, xName, yName, std::string("Channel"));
        return;
    }

    if (flag == LinacResults[4])
    {
        yName = "X, m";
        xName = "dX/dZ, rad";

        std::vector<std::vector<float>> Xar(numberOfTraces);
        std::vector<std::vector<float>> Yar1(numberOfTraces);
        std::vector<std::vector<float>> Yar2(numberOfTraces);

        for (int k = 0; k < numberOfTraces; k++)
        {

            int    tr = rand() % data->data[2].size();
            double V  = data->data[1][0].back();

            double em = data->emittanceX * 1e-6;

            double betaX  = data->data[2][tr].back() / (em);
            double AlphaX = -data->data[3][tr].back() / em;
            double gammaX = data->data[4][tr].back() / (em);

            double xM  = sqrt(data->data[2][tr].back());
            double dxM = sqrt(data->data[4][tr].back());
            int    n   = 500;
            double dx  = 2 * xM / n;
            double xc  = -1.5 * xM;
            double dxC;
            while (1)
            {
                if (std::abs(xM - std::abs(xc)) < 0.1 * xM)
                    dxC = 0.04 * dx;
                else
                    dxC = dx;

                Xar[k].push_back(xc);

                double b = 2 * AlphaX * xc;
                double d = b * b - 4 * betaX * (gammaX * xc * xc - em);

                Yar1[k].push_back((1.0 / V) * (-b + sqrt(d)) / (2 * betaX));
                Yar2[k].push_back((1.0 / V) * (-b - sqrt(d)) / (2 * betaX));

                xc = xc + dxC;

                if (xc > 1.5 * xM)
                    break;
            }
        }
        input->addVisualizationDataPlotsAutoY(Xar, Yar1, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(Xar, Yar2, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        return;
    }

    if (flag == LinacResults[5])
    {
        yName = "Y, m";
        xName = "dY/dZ, rad";

        std::vector<std::vector<float>> Xar(numberOfTraces);
        std::vector<std::vector<float>> Yar1(numberOfTraces);
        std::vector<std::vector<float>> Yar2(numberOfTraces);

        for (int k = 0; k < numberOfTraces; k++)
        {

            int    tr = rand() % data->data[5].size();
            double V  = data->data[1][0].back();

            double em = data->emittanceX * 1e-6;

            double betaX  = data->data[5][tr].back() / (em);
            double AlphaX = -data->data[6][tr].back() / em;
            double gammaX = data->data[7][tr].back() / (em);

            double xM  = sqrt(data->data[5][tr].back());
            double dxM = sqrt(data->data[7][tr].back());
            int    n   = 500;
            double dx  = 2 * xM / n;
            double xc  = -1.5 * xM;
            double dxC;
            while (1)
            {
                if (std::abs(xM - std::abs(xc)) < 0.1 * xM)
                    dxC = 0.04 * dx;
                else
                    dxC = dx;

                Xar[k].push_back(xc);

                double b = 2 * AlphaX * xc;
                double d = b * b - 4 * betaX * (gammaX * xc * xc - em);

                Yar1[k].push_back((1.0 / V) * (-b + sqrt(d)) / (2 * betaX));
                Yar2[k].push_back((1.0 / V) * (-b - sqrt(d)) / (2 * betaX));

                xc = xc + dxC;

                if (xc > 1.5 * xM)
                    break;
            }
        }
        input->addVisualizationDataPlotsAutoY(Xar, Yar1, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(Xar, Yar2, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        return;
    }

    if (flag == LinacResults[6])
    {
        yName = "Transmission";
        xName = "T, ns";
        ymax  = 1.0;
        ymin  = 1.0;
        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[2], 2, xName, yName,
                                        std::string("Transmission"), ymax, ymin);
        yName = "Acceleration";
        xName = "T, ns";
        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[3], 2, xName, yName,
                                        std::string("Acceleration"), ymax, ymin);
        return;
    }

    /*
    yName = "R chan";	xName = "T, ns";
    input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[0], 2, xName, yName, std::string("R channel"),
    ymax, ymin); yName = "R beam";	xName = "T, ns"; input->addVisualizationDataPlot(data->TimeArrayAdd,
    data->dataAdd[1], 2, xName, yName, std::string("R beam"), ymax, ymin); return;

    case 8:
                    yName = "Transmission";	xName = "T, ns";
                    input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[2], 2, xName, yName,
    std::string("Transmission"), ymax, ymin); yName = "Acceleration";	xName = "T, ns";
                    input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[3], 2, xName, yName,
    std::string("Acceleration"), ymax, ymin); break;

            case 9:
                    yName = "Defocusing";	xName = "T, ns";
                    input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[5], 2, xName, yName,
    std::string("Defocusing"), ymax, ymin); break;


            case 10:
                    yName = "Energy, Ev";	xName = "Z, m";
                    std::vector<float> tmp = data->data[3][0];
                    std::vector<float> tmpT = data->TimeArray[0];
                    std::vector<float> tmpT1 = data->dataAdd[14];

                    en = LIGHT_VELOCITY*LIGHT_VELOCITY*data->mass;

                    for (int i = 0; i < tmp.size(); i++)
                            tmp[i] = tmp[i] * LIGHT_VELOCITY;

                    float omega = 2 * PI*108e6;

                    for (int i = 0; i < tmp.size(); i++)
                    {
                            tmpT[i] = tmpT[i] * 1e-9*omega;
                            while (tmpT[i]> PI)
                                    tmpT[i] = tmpT[i] - 2 * PI;
                    }

                    std::vector<float> cells;
                    for (int i = 0; i < data->dataAdd[14].size(); i++)
                    {
                            cells.push_back(i);

                            while (tmpT1[i]> PI / 2)
                                    tmpT1[i] = tmpT1[i] - PI;
                    }
            input->addVisualizationDataPlot(cells, tmpT1, 2, xName, yName, std::string("Velocity"), ymax, ymin);

            }*/
};
