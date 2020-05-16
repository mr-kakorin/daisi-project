#include <algorithm>
#include <memory>
#include <numeric>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>
#include <daisi-solver/FlagStringsSolver.h>
#include <daisi-solver/SynchrotronStrings.h>

#include "GeneralTools.h"
#include "VisualizationFunctionsNucl.h"
#include "colors.h"
#include "ellipse.h"
#include "vtkComponent.h"

void GetParticlesCloud(std::vector<double> plotProperties, const std::vector<std::vector<float>>& TimeArray,
                       const std::vector<std::vector<float>>& x1, std::vector<std::vector<float>>& x2,
                       std::vector<float>& x, std::vector<float>& y)
{
    int j             = 0;
    int iplot         = TimeArray[0].size() - 1;
    int TimeArrayplot = 0;
    for (int k = 0; k < TimeArray.size(); k++)
    {
        if (TimeArray[k].size() - 1 > iplot)
        {
            TimeArrayplot = k;
            iplot         = TimeArray[k].size() - 1;
        }
    }

    if (plotProperties[1] < 0)
    {
        iplot = TimeArray[TimeArrayplot].size() - 1;
    }
    else
    {
        for (j = 0; j < TimeArray[TimeArrayplot].size() - 1; j++)
        {
            if (plotProperties[1] >= TimeArray[0][j] && plotProperties[1] <= TimeArray[0][j + 1])
                break;
        };
        if (std::abs(TimeArray[0][j] - plotProperties[1]) < std::abs(TimeArray[0][j + 1] - plotProperties[1]))
            iplot = j;
        else
            iplot = j + 1;
    }

    int ii = 0;

    for (int i = 0; i < x1.size(); i++)
    {
        if (iplot < x1[i].size())
        {
            x[ii] = x1[i][iplot];
            y[ii] = x2[i][iplot];
            ii++;
        }
        else
        {
            x.pop_back();
            y.pop_back();
        };
    };
}
void ShowSimulationNuclTwiss(vtkComponent* input, AccelModelInterface* currentModel,
                             std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                             std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                             float& ymin)
{

    int maxPoint = simData->XData[0].size();
    if (plotProperties[1] >= 0)
    {
        for (int i = 0; i < maxPoint; i++)
        {
            if (simData->XData[0][i] > plotProperties[1])
            {
                maxPoint = i;
                break;
            };
        }
    };

    std::string yName;
    std::string xName = "L, m";

    if (flag == SynchrotronTwissDynFlags[0])
    {
        yName = "BETA";
        input->addVisualizationDataPlot(simData->XData[0], simData->YData[0][0], 2.0, xName, yName, "BETAx", ymax, ymin,
                                        false, 0, maxPoint);
        input->addVisualizationDataPlot(simData->XData[3], simData->YData[3][0], 2.0, xName, yName, "BETAy", ymax, ymin,
                                        false, 0, maxPoint);
    }
    if (flag == SynchrotronTwissDynFlags[1])
    {
        yName = "ALFA";
        input->addVisualizationDataPlot(simData->XData[1], simData->YData[1][0], 2.0, xName, yName, "ALFAx", ymax, ymin,
                                        false, 0, maxPoint);
        input->addVisualizationDataPlot(simData->XData[4], simData->YData[4][0], 2.0, xName, yName, "ALFAy", ymax, ymin,
                                        false, 0, maxPoint);
    }
    if (flag == SynchrotronTwissDynFlags[2])
    {
        yName = "MU";
        input->addVisualizationDataPlot(simData->XData[2], simData->YData[2][0], 2.0, xName, yName, "MUx", ymax, ymin,
                                        false, 0, maxPoint);
        input->addVisualizationDataPlot(simData->XData[5], simData->YData[5][0], 2.0, xName, yName, "MUy", ymax, ymin,
                                        false, 0, maxPoint);
    }
};
void ShowSimulationNuclBeam(vtkComponent* input, AccelModelInterface* currentModel,
                            std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                            std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                            float& ymin)
{

    std::string yName;
    std::string xName = "L, m";

    if (!simData->TimeArray.size())
        return;

    int maxPoint = simData->TimeArray[0].size();
    if (plotProperties[1] >= 0)
    {
        for (int i = 0; i < maxPoint; i++)
        {
            if (simData->TimeArray[0][i] > plotProperties[1])
            {
                maxPoint = i;
                break;
            };
        }
    };

    if (flag == SynchrotronDynFlags[0])
    {
        yName = "X, m";

        input->addVisualizationDataPlotsAutoY(simData->TimeArray, simData->data[0], colors::blackcolor, 1, xName, yName,
                                              plotProperties[0], 0, props, ymax, ymin, 0, maxPoint);
        if (simData->TimeArrayAdd.size())
            input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[0], colors::redcolor, 2, "", "", 0,
                                            1, 0, "Aperture", 1, false);
    }
    if (flag == SynchrotronDynFlags[1])
    {
        yName = "dX/dZ, rad";
        input->addVisualizationDataPlotsAutoY(simData->TimeArray, simData->data[1], colors::blackcolor, 1, xName, yName,
                                              plotProperties[0], 0, props, ymax, ymin, 0, maxPoint);
    }
    if (flag == SynchrotronDynFlags[2])
    {
        yName = "Y, m";
        input->addVisualizationDataPlotsAutoY(simData->TimeArray, simData->data[2], colors::blackcolor, 1, xName, yName,
                                              plotProperties[0], 0, props, ymax, ymin, 0, maxPoint);
        if (simData->TimeArrayAdd.size())
            input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[1], colors::redcolor, 2, "", "", 0,
                                            1, 0, "Aperture", 1, false);
    }
    if (flag == SynchrotronDynFlags[3])
    {
        yName = "dY/dZ, rad";
        input->addVisualizationDataPlotsAutoY(simData->TimeArray, simData->data[3], colors::blackcolor, 1, xName, yName,
                                              plotProperties[0], 0, props, ymax, ymin, 0, maxPoint);
    }
    if (flag == SynchrotronDynFlags[4])
    {
        xName   = "X, m";
        yName   = "dX/dZ, rad";
        flagSeq = 1;
        std::vector<float> x(simData->data[2].size());
        std::vector<float> y(simData->data[2].size());

        GetParticlesCloud(plotProperties, simData->TimeArray, simData->data[0], simData->data[1], x, y);

        input->addVisualizationDataPlot(x, y, colors::blackcolor, 4.0, xName, yName, 1, 1, 0, "XdX plane distribution",
                                        1.0, false);

        std::vector<std::vector<float>> x1El;
        std::vector<std::vector<float>> x2El;

        std::vector<std::vector<float>> yEl1;
        std::vector<std::vector<float>> yEl2;
        std::vector<std::string>        out;

        calculateEllipsesParameters(out, x1El, x2El, yEl1, yEl2, x, y, {float(0.99), float(0.9)});

        for (int i = 0; i < out.size(); i++)
        {
            input->addVisualizationDataPlot(x1El[i], yEl1[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
            input->addVisualizationDataPlot(x2El[i], yEl2[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
        }

        input->AddText(out, colors::colors, 250);
        return;
    }

    if (flag == SynchrotronDynFlags[5])
    {
        int j   = 0;
        xName   = "Y, m";
        yName   = "dY/dZ, rad";
        flagSeq = 1;
        std::vector<float> x(simData->data[2].size());
        std::vector<float> y(simData->data[2].size());

        GetParticlesCloud(plotProperties, simData->TimeArray, simData->data[2], simData->data[3], x, y);

        input->addVisualizationDataPlot(x, y, colors::blackcolor, 4.0, xName, yName, 1, 1, 0, "YdY plane distribution",
                                        1.0, false);

        std::vector<std::vector<float>> x1El;
        std::vector<std::vector<float>> x2El;

        std::vector<std::vector<float>> yEl1;
        std::vector<std::vector<float>> yEl2;
        std::vector<std::string>        out;

        calculateEllipsesParameters(out, x1El, x2El, yEl1, yEl2, x, y, {float(0.99), float(0.9)});

        for (int i = 0; i < out.size(); i++)
        {
            input->addVisualizationDataPlot(x1El[i], yEl1[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
            input->addVisualizationDataPlot(x2El[i], yEl2[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
        }

        input->AddText(out, colors::colors, 250);

        return;
    }
    if (flag == SynchrotronDynFlags[6])
    {
        int j   = 0;
        xName   = "X, m";
        yName   = "Y, m";
        flagSeq = 1;
        std::vector<float> x(simData->data[2].size());
        std::vector<float> y(simData->data[2].size());

        GetParticlesCloud(plotProperties, simData->TimeArray, simData->data[0], simData->data[2], x, y);

        input->addVisualizationDataPlot(x, y, colors::blackcolor, 4.0, xName, yName, 1, 1, 0, "XY plane distribution",
                                        1.0, false);

        return;
    }
    if (flag == SynchrotronDynFlags[7])
    {
        flagSeq = 1;
        if (simData->TimeArrayAdd.size())
            input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[2], colors::redcolor, 2, "L, m",
                                            "Transmission", 0, 1, 0, "Transmission", 1, false);
        return;
    }
    if (flag == SynchrotronDynFlags[8])
    {
        float target = 0.9;
        if (plotProperties[0] > 0.01 && plotProperties[0] < 1.0)
            target = plotProperties[0];
        auto alpha_bet_x = calc_alpha_beta(simData->data[0], simData->data[1], target);
        auto alpha_bet_y = calc_alpha_beta(simData->data[2], simData->data[3], target);

        yName = "BETA";
        input->addVisualizationDataPlot(simData->TimeArray[0], alpha_bet_x.second, 2.0, xName, yName, "BETAx", ymax,
                                        ymin, false, 0, maxPoint);
        input->addVisualizationDataPlot(simData->TimeArray[0], alpha_bet_y.second, 2.0, xName, yName, "BETAy", ymax,
                                        ymin, false, 0, maxPoint);
    }
    if (flag == SynchrotronDynFlags[9])
    {
        float target = 0.9;
        if (plotProperties[0] > 0.01 && plotProperties[0] < 1.0)
            target = plotProperties[0];
        auto alpha_bet_x = calc_alpha_beta(simData->data[0], simData->data[1], target);
        auto alpha_bet_y = calc_alpha_beta(simData->data[2], simData->data[3], target);

        yName = "Alpha";
        input->addVisualizationDataPlot(simData->TimeArray[0], alpha_bet_x.first, 2.0, xName, yName, "Alphax", ymax,
                                        ymin, false, 0, maxPoint);
        input->addVisualizationDataPlot(simData->TimeArray[0], alpha_bet_y.first, 2.0, xName, yName, "Alphay", ymax,
                                        ymin, false, 0, maxPoint);
    }
    if (flag == SynchrotronDynFlags[10])
    {
        auto cm_x = calc_cm(simData->data[0], simData->data[1]);
        auto cm_y = calc_cm(simData->data[2], simData->data[3]);

        yName = "Center of mass";

        input->addVisualizationDataPlot(simData->TimeArray[0], cm_x.first, 2.0, xName, yName, "X", ymax, ymin, false, 0,
                                        maxPoint);
        input->addVisualizationDataPlot(simData->TimeArray[0], cm_y.first, 2.0, xName, yName, "Y", ymax, ymin, false, 0,
                                        maxPoint);
    }
    if (flag == SynchrotronDynFlags[11])
    {
        auto cm_x = calc_cm(simData->data[0], simData->data[1]);
        auto cm_y = calc_cm(simData->data[2], simData->data[3]);

        yName = "Center of mass";

        input->addVisualizationDataPlot(simData->TimeArray[0], cm_x.second, 2.0, xName, yName, "X", ymax, ymin, false,
                                        0, maxPoint);
        input->addVisualizationDataPlot(simData->TimeArray[0], cm_y.second, 2.0, xName, yName, "Y", ymax, ymin, false,
                                        0, maxPoint);
    }
};

void ShowNuclElements(vtkComponent* input, AccelModelInterface* currentModel,
                      std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                      std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin)
{
    float* color = NULL;
    color        = new float[3];

    std::vector<std::string> names = simData->names;

    std::vector<std::vector<double>> propSeqs = simData->props;

    if (!propSeqs.size())
        return;

    int maxPoint = propSeqs[0].size();
    if (plotProperties[1] >= 0)
    {
        for (int i = 0; i < maxPoint; i++)
        {
            if (propSeqs[0][i] > plotProperties[1])
            {
                maxPoint = i;
                break;
            };
        }
    };

    float ys = ymax + (ymax - ymin) / 4;
    float ye = ys + (ymax - ymin) / 4;

    for (int i = 0; i < maxPoint; i++)
    {

        color[0] = -1;
        //"RBEND", "QUADRUPOLE", "KICKER"
        if (names[i] == "RBEND")
        {
            color[0] = colors::bluecolor[0];
            color[1] = colors::bluecolor[1];
            color[2] = colors::bluecolor[2];
        }
        if (names[i] == "KICKER")
        {
            color[0] = colors::greencolor[0];
            color[1] = colors::greencolor[1];
            color[2] = colors::greencolor[2];
        }
        if (names[i] == "QUADRUPOLE")
        {
            color[0] = colors::redcolor[0];
            color[1] = colors::redcolor[1];
            color[2] = colors::redcolor[2];
        }
        if (color[0] != -1)
        {
            input->addVisualizationDataPlot({float(propSeqs[0][i]), float(propSeqs[0][i] + propSeqs[1][i])}, {ys, +ys},
                                            color, 2, "", "", 0, 1, 0, names[i], 1, false);
            input->addVisualizationDataPlot({float(propSeqs[0][i]), float(propSeqs[0][i] + propSeqs[1][i])}, {ye, ye},
                                            color, 2, "", "", 0, 1, 0, names[i], 1, false);
            input->addVisualizationDataPlot({float(propSeqs[0][i]), float(propSeqs[0][i])}, {ys, ye}, color, 2, "", "",
                                            0, 1, 0, names[i], 1, false);
            input->addVisualizationDataPlot(
                {float(propSeqs[0][i] + propSeqs[1][i]), float(propSeqs[0][i] + propSeqs[1][i])}, {ys, ye}, color, 2,
                "", "", 0, 1, 0, names[i], 1, false);
            input->addVisualizationDataPlot({float(propSeqs[0][i]), float(propSeqs[0][i] + propSeqs[1][i])}, {ys, ye},
                                            color, 2, "", "", 0, 2, 0, names[i], 1, false);
            input->addVisualizationDataPlot({float(propSeqs[0][i] + propSeqs[1][i]), float(propSeqs[0][i])}, {ys, ye},
                                            color, 2, "", "", 0, 2, 0, names[i], 1, false);
        }
        // input->addVisualizationDataPlot( { float(propSeqs[0][i]),float(propSeqs[0][i] + propSeqs[1][i]) }, { 0.0, 0.0
        // }, color, 140, "", "", 0, 1, 0, names[i]);
    }
};

void ShowNuclCorrOrbit(vtkComponent* input, AccelModelInterface* currentModel,
                       std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                       std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin)
{
    std::string yName;
    std::string xName = "Pick-up electrode number";

    if (flag == SynchrotronCorrFlags[0])
    {
        std::vector<float> numbers;
        for (int i = 0; i < simData->YData[4][0].size(); i++)
            numbers.push_back(i);

        yName = "I, rad";
        xName = "Corrector number";
        input->addVisualizationDataPlot(numbers, simData->YData[4][0], 2.0, xName, yName,
                                        "XdX plane correction current", ymax, ymin, false, 0);
        input->addVisualizationDataPlot(numbers, simData->YData[5][0], 2.0, xName, yName,
                                        "YdY plane correction current", ymax, ymin, false, 0);
    }
    if (flag == SynchrotronCorrFlags[1])
    {
        std::vector<float> numbers;
        for (int i = 0; i < simData->YData[0][0].size(); i++)
            numbers.push_back(i);

        yName = "X coordinate deviatian";
        input->addVisualizationDataPlot(numbers, simData->YData[0][0], 2.0, xName, yName, "Deviation after correction",
                                        ymax, ymin, false, 0);
        input->addVisualizationDataPlot(numbers, simData->YData[2][0], 2.0, xName, yName, "Initial deviation", ymax,
                                        ymin, false, 0);
    }
    if (flag == SynchrotronCorrFlags[2])
    {
        std::vector<float> numbers;
        for (int i = 0; i < simData->YData[1][0].size(); i++)
            numbers.push_back(i);

        yName = "Y coordinate deviatian";
        input->addVisualizationDataPlot(numbers, simData->YData[1][0], 2.0, xName, yName, "Deviation after correction",
                                        ymax, ymin, false, 0);
        input->addVisualizationDataPlot(numbers, simData->YData[3][0], 2.0, xName, yName, "Initial deviation", ymax,
                                        ymin, false, 0);
    }
    if (flag == SynchrotronCorrFlags[3])
    {
        std::vector<float> numbers;
        for (int i = 0; i < simData->YData[1][0].size(); i++)
            numbers.push_back(i);

        yName = "Singular value";
        xName = "Number of singular value";

        input->addVisualizationDataPlot(numbers, simData->YData[6][0], 2.0, xName, yName, "Singular value for X matrix",
                                        ymax, ymin, false, 0);
        input->addVisualizationDataPlot(numbers, simData->YData[7][0], 2.0, xName, yName, "Singular value for Y matrix",
                                        ymax, ymin, false, 0);
    }
};

void ShowSimulationNuclCM(vtkComponent* input, AccelModelInterface* currentModel,
                          std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                          std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin)
{
    std::string xName = "L, m";
    std::string yName;

    if (flag == SynchrotronDynCMFlags[0])
    {
        yName = "X and Y, m";
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[0][1], 2.0, xName, yName, "X", ymax, ymin,
                                        false);
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[0][3], 2.0, xName, yName, "Y", ymax, ymin,
                                        false);

        if (simData->TimeArrayAdd.size())
        {
            input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[0], colors::blackcolor, 2, "", "",
                                            0, 1, 0, "Aperture X", 1, false);
            input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[1], colors::redcolor, 2, "", "", 0,
                                            1, 0, "Aperture Y", 1, false);
            ymax = simData->dataAdd[0][0];
        }
    }
    if (flag == SynchrotronDynCMFlags[1])
    {
        yName = "dX/dZ and dY/dZ, rad";
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[0][2], 2.0, xName, yName, "dX/dZ", ymax,
                                        ymin, false);
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[0][4], 2.0, xName, yName, "dY/dZ", ymax,
                                        ymin, false);
    }
};
void ShowSimulationNuclOrbit(vtkComponent* input, AccelModelInterface* currentModel,
                             std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                             std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                             float& ymin)
{
    std::string xName = "L, m";
    std::string yName;

    if (flag == SynchrotronOrbitFlags[0])
    {
        flagSeq = 0;
        yName   = "X and Y, m";
        input->addVisualizationDataPlot(simData->XData[0], simData->YData[0][0], 2.0, xName, yName, "X", ymax, ymin,
                                        false);
        input->addVisualizationDataPlot(simData->XData[1], simData->YData[1][0], 2.0, xName, yName, "Y", ymax, ymin,
                                        false);
        input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[0], colors::blackcolor, 2, "", "", 0, 1,
                                        0, "Aperture X", 1, false);
        input->addVisualizationDataPlot(simData->TimeArrayAdd, simData->dataAdd[1], colors::redcolor, 2, "", "", 0, 1,
                                        0, "Aperture Y", 1, false);
        ymax = simData->dataAdd[0][0];
    }
};

void ShowSimulationNuclAcceptance(vtkComponent* input, AccelModelInterface* currentModel,
                                  std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                                  std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                                  float& ymin)
{
    std::string yName;
    std::string xName = "L, m";
    if (flag == SynchrotronAccFlags[0])
    {
        xName   = "X, m";
        yName   = "dX/dZ, rad";
        flagSeq = 1;
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[1][0], colors::blackcolor, 4.0, xName,
                                        yName, 1, 1, 0, "XdX plane distribution", 1.0, false);

        std::vector<std::vector<float>> x1El;
        std::vector<std::vector<float>> x2El;

        std::vector<std::vector<float>> yEl1;
        std::vector<std::vector<float>> yEl2;
        std::vector<std::string>        out;

        calculateEllipsesParameters(out, x1El, x2El, yEl1, yEl2, simData->YData[0][0], simData->YData[1][0],
                                    {float(0.97)});

        for (int i = 0; i < out.size(); i++)
        {
            input->addVisualizationDataPlot(x1El[i], yEl1[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
            input->addVisualizationDataPlot(x2El[i], yEl2[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
        }

        input->AddText(out, colors::colors, 250);
        return;
    }
    if (flag == SynchrotronAccFlags[1])
    {
        xName   = "Y, m";
        yName   = "dY/dZ, rad";
        flagSeq = 1;
        input->addVisualizationDataPlot(simData->YData[2][0], simData->YData[3][0], colors::blackcolor, 4.0, xName,
                                        yName, 1, 1, 0, "XdX plane distribution", 1.0, false);

        std::vector<std::vector<float>> x1El;
        std::vector<std::vector<float>> x2El;

        std::vector<std::vector<float>> yEl1;
        std::vector<std::vector<float>> yEl2;
        std::vector<std::string>        out;

        calculateEllipsesParameters(out, x1El, x2El, yEl1, yEl2, simData->YData[2][0], simData->YData[3][0],
                                    {float(0.97)});

        for (int i = 0; i < out.size(); i++)
        {
            input->addVisualizationDataPlot(x1El[i], yEl1[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
            input->addVisualizationDataPlot(x2El[i], yEl2[i], colors::colors[i], 2, xName, yName, 0, 1, 0, "Ellipse",
                                            1.0, false);
        }

        input->AddText(out, colors::colors, 250);
        return;
    }

    if (flag == SynchrotronAccFlags[2])
    {
        xName   = "X, m";
        yName   = "Y, m";
        flagSeq = 1;
        input->addVisualizationDataPlot(simData->YData[0][0], simData->YData[2][0], colors::blackcolor, 4.0, xName,
                                        yName, 1, 1, 0, "XdX plane distribution", 1.0, false);

        return;
    }
};

void SynchrotronRec(vtkComponent* input, AccelModelInterface* currentModel,
                    std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                    std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin)
{
    std::string yName;
    std::string xName = "L, m";

    auto lambda = [&](int indShow, std::string yName, int offSet) {

        float fit = 0;
        float min = 1e38;
        int   indMin;
        for (size_t i = 4 + offSet; i <= simData->YData.size() - 1; i = i + 2)
        {
            fit = 0;
            for (size_t j = 0; j != simData->YData[i][indShow].size(); j++)
            {
                float sq = (simData->YData[2 + offSet][indShow][j] - simData->YData[i][indShow][j]);
                fit      = fit + sq;
            }
            if (fit < min)
            {
                min    = fit;
                indMin = i;
            }
        }

        input->addVisualizationDataPlot(simData->YData[0 + offSet][0], simData->YData[0 + offSet][indShow],
                                        colors::blackcolor, 4, xName, yName, 0, 1, 0, "hgfhg", 1, false);

        input->addVisualizationDataPlot(simData->YData[2 + offSet][0], simData->YData[2 + offSet][indShow],
                                        colors::greencolor, 4, xName, yName, 0, 1, 0, "hgfhg", 1, false);
        for (size_t ind = 4 + offSet; ind <= simData->YData.size() - 1; ind = ind + 2)
        {
            auto color = colors::redcolor;
            if (indMin == ind)
                color = colors::bluecolor;

            input->addVisualizationDataPlot(simData->YData[ind][0], simData->YData[ind][indShow], color, 2, xName,
                                            yName, 0, 1, 0, "hgfhg", 1, false);
        }
        std::vector<std::string> out;
        input->AddText(
            {"Non-dev. sol.", "Dev. sol.", "Rec. sol.",
             "Best rec. sol. " + std::to_string(simData->XData.back()[(indMin - 4 - offSet) / 2]) + " monitors"},
            {colors::blackcolor, colors::greencolor, colors::redcolor, colors::bluecolor}, 100);
    };

    if (flag == SynchrotronRecCMFlags[0])
    {
        lambda(1, SynchrotronRecCMFlags[0], 0);
    }
    if (flag == SynchrotronRecCMFlags[1])
    {
        lambda(2, SynchrotronRecCMFlags[1], 0);
    }
    if (flag == SynchrotronRecCMFlags[2])
    {
        lambda(3, SynchrotronRecCMFlags[2], 0);
    }
    if (flag == SynchrotronRecCMFlags[3])
    {
        lambda(4, SynchrotronRecCMFlags[3], 0);
    }
    if (flag == SynchrotronRecCMFlags[4])
    {
        for (size_t ind = 0; ind != simData->XData.size(); ind++)
        {

            std::vector<float> x;

            for (int s = 0; s < simData->XData[ind].size(); s++)
                x.push_back(s);

            input->addVisualizationDataPlot(x, simData->XData[ind], colors::blackcolor, 2, "Iteration number",
                                            "Fitness", 0, 1, 0, "Aperture X", 1, false);
        }
    }
    if (flag == SynchrotronRecCMFlags[5])
    {
        lambda(1, SynchrotronRecCMFlags[5], 1);
    }
    if (flag == SynchrotronRecCMFlags[6])
    {
        lambda(2, SynchrotronRecCMFlags[6], 1);
    }
    if (flag == SynchrotronRecCMFlags[7])
    {
        lambda(4, SynchrotronRecCMFlags[7], 1);
    }
    if (flag == SynchrotronRecCMFlags[8])
    {
        lambda(5, SynchrotronRecCMFlags[8], 1);
    }
    if (flag == SynchrotronRecCMFlags[9])
    {
        std::vector<float> x;
        std::vector<float> y;

        if (plotProperties[0] > simData->YData1[0].size() - 1)
        {
            return;
        }

        for (int s = 0; s < simData->YData1[0][plotProperties[0]].size(); s++)
        {
            x.push_back(s);
            y.push_back(simData->YData1[1][0][plotProperties[0]]);
        }

        input->addVisualizationDataPlot(x, simData->YData1[0][plotProperties[0]], colors::blackcolor, 2,
                                        "Iteration number", "Parameter", 0, 1, 0, "sads", 1, false);

        input->addVisualizationDataPlot(x, y, colors::redcolor, 2, "Iteration number", "Parameter finish", 0, 1, 0,
                                        "sads", 1, false);
    }
}
void SynchrotronShuffle(vtkComponent* input, AccelModelInterface* currentModel,
	std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
	std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
	float& ymin)
{
	std::vector<float> XChart;
	std::vector<float> YChart;   
	std::vector<float> XChart_line;
	std::vector<float> YChart_line;
	std::vector<float> YChart_gauss;
	std::vector<float> XChart_gauss;
	int    nCharts = 30;

	auto plot = [&](const int p) {
		if (!getChart(simData->XData[p], nCharts, XChart, YChart, 0, -1))
		{
			return;
		}

		YChart_gauss = YChart;
		XChart_gauss = XChart;
		float dx_ch = XChart[1] - XChart[0];
		std::for_each(YChart_gauss.begin(), YChart_gauss.end(), [&dx_ch](auto& y) { y = y / (100.0 * dx_ch); });
		std::for_each(XChart_gauss.begin(), XChart_gauss.end(), [&dx_ch](auto& y) { y = y + dx_ch; });

		convertChart(XChart, YChart, XChart_line, YChart_line);

		auto gauss_par = estimate_gauss(XChart_gauss, YChart_gauss, 1e-3);

		convertChart(XChart, YChart, XChart_line, YChart_line);
		input->addVisualizationDataPlot(XChart_line, YChart_line, 4.0, "X_max", "% of experiments", "% of experiments", ymax, ymin,
			false);

		input->AddText({ "mu = " + std::to_string(gauss_par[0]), "sigma = " + std::to_string(gauss_par[1]) },
		{ colors::blackcolor, colors::blackcolor }, 100);
	};

	if (flag == SynchrotronShuffleFlags[0])
	{
		plot(0);
	}
	if (flag == SynchrotronShuffleFlags[1])
	{
		plot(1);
	}
}

void SynchrotronEst(vtkComponent* input, AccelModelInterface* currentModel,
                    std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                    std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin)
{

    auto lambda = [&](int N, std::string xName, std::string yName, int flag) {

        size_t n_div = 0;
        if (plotProperties[1] >= 0)
        {
            n_div = plotProperties[1];
        }
        double wmax_in = -1;
        int    nCharts = 30;

        std::vector<float> XChart;
        std::vector<float> YChart;
        std::vector<float> YChart_gauss;
        std::vector<float> XChart_gauss;

        std::vector<float> XChart_line;
        std::vector<float> YChart_line;
        std::vector<float> gauss_par;

        std::vector<float> xval = simData->XData[N];

        if (N == 2 || N == 3 || N == 6 || N == 7)
        {
            xval.clear();
            int NN = N + 2;
            if (N == 2 || N == 3)
            {
                NN = N + 6;
            }
            std::vector<size_t> indexes(simData->XData[N].size());
            std::iota(indexes.begin(), indexes.end(), 0);
            float av = std::accumulate(simData->XData[NN].begin(), simData->XData[NN].end(), 0.0) /
                       float(simData->XData[8].size());

            for (size_t div = 0; div < n_div; div++)
            {
                indexes.clear();
                for (size_t i = 0; i < simData->XData[NN].size(); i++)
                {
                    if (simData->XData[NN][i] < av)
                    {
                        indexes.push_back(i);
                    }
                }
                av = 0;
                for (size_t i = 0; i < indexes.size(); i++)
                {
                    av = av + simData->XData[NN][indexes[i]];
                }
                av = av / indexes.size();
            }
            for (size_t i = 0; i < indexes.size(); i++)
            {
                xval.push_back(simData->XData[N][indexes[i]]);
            }
            n_div = 0;
        }

        for (size_t div = 0; div < n_div + 1; div++)
        {
            YChart.clear();
            XChart.clear();
            YChart_gauss.clear();
            XChart_gauss.clear();
            XChart_line.clear();
            YChart_line.clear();
            gauss_par.clear();

            if (!getChart(xval, nCharts, XChart, YChart, flag, wmax_in))
            {
                return;
            }

            YChart_gauss = YChart;
            XChart_gauss = XChart;
            float dx_ch  = XChart[1] - XChart[0];
            std::for_each(YChart_gauss.begin(), YChart_gauss.end(), [&dx_ch](auto& y) { y = y / (100.0 * dx_ch); });
            std::for_each(XChart_gauss.begin(), XChart_gauss.end(), [&dx_ch](auto& y) { y = y + dx_ch; });

            convertChart(XChart, YChart, XChart_line, YChart_line);

            gauss_par = estimate_gauss(XChart_gauss, YChart_gauss, 1e-3);

            wmax_in = gauss_par[0];
        }

        input->addVisualizationDataPlot(XChart_line, YChart_line, 4.0, xName, yName, "% of experiments", ymax, ymin,
                                        false);

        input->AddText({"mu = " + std::to_string(gauss_par[0]), "sigma = " + std::to_string(gauss_par[1])},
                       {colors::blackcolor, colors::blackcolor}, 100);
        // input->addVisualizationDataPlot(XChart, YChart, colors::blackcolor, 5, xName, yName, 2, 1, 0, "distribution",
        //	1.0, false);
        /*float yMin; float yMax;
        input->addVisualizationDataPlots(XChart_line,
                YChart_line, colors::blackcolor,
                1.0, xName, yName, 0, 0,
                std::vector<int>{}, yMin, yMax);*/
    };

    size_t k = 0;
    for (const auto& s : simData->dataFlags)
    {
        if (s == flag)
        {
            if (k == simData->dataFlags.size() - 2)
            {
                lambda(6, flag, "% of experiments", 1);
                break;
            }
            if (k == simData->dataFlags.size() - 1)
            {
                lambda(7, flag, "% of experiments", 1);
                break;
            }
            lambda(k, flag, "% of experiments", 0);
            break;
        }
        k++;
    }
    /*if (flag == SynchrotronEstFlags[0])
    {
            lambda(0, "BetaX max", "Number of experiments");
    }
    if (flag == SynchrotronEstFlags[1])
    {
            lambda(1, "BetaY max", "Number of experiments");

    }*/
}