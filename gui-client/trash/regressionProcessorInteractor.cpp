#include "DaisiClient.h"
#include "VisualizationFunctions.h"
#include "regressionProcessor.h"
#include "vtkComponent.h"
void Daizy::regressionProcessor()
{

    progress  = 0;
    flagAbort = true;

    std::vector<std::vector<float>> resultsLineX;
    std::vector<std::vector<float>> resultsLineY;

    std::vector<std::vector<float>> results;
    std::vector<std::vector<float>> eps;

    processOilwells(errorMsg, progress, results, resultsLineX, resultsLineY, eps);
    AddResultPlot();
    AddResultPlot();
    AddResultPlot();

    VTKArray[0]->setDataFunction(SimpleShowPoints, results[0], results[1],
                                 std::vector<std::vector<float>>{resultsLineX.begin(), resultsLineX.begin() + 2},
                                 std::vector<std::vector<float>>{resultsLineY.begin(), resultsLineY.begin() + 2});
    VTKArray[0]->refresh(0);

    VTKArray[1]->setDataFunction(SimpleShowPoints, results[0], results[2],
                                 std::vector<std::vector<float>>{resultsLineX.begin() + 2, resultsLineX.end()},
                                 std::vector<std::vector<float>>{resultsLineY.begin() + 2, resultsLineY.end()});

    VTKArray[1]->refresh(0);

    VTKArray[2]->setDataFunction(SimpleShowChart, eps[0]);

    VTKArray[2]->refresh(0);

    VTKArray[3]->setDataFunction(SimpleShowChart, eps[1]);

    VTKArray[3]->refresh(0);
};