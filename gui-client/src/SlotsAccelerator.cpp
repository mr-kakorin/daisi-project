#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"

void Daizy::AccelSolve()
{
    // connect(timerViz, SIGNAL(timeout()), this, SLOT(updateGrViz()));

    // timerViz->start(1000);

    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;

    bool ok = middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    if (!ok)
    {
        QMessageBox::critical(this, "Daisi error", "Non-digit input is converted to zero");
        return;
    }

    std::string         errorMessage;
    std::vector<double> parameters22;
    for (int i = 0; i < parameters1.size() - 1; i++)
        parameters22.push_back(parameters1[1 + i][0]);

    currentProject->accelModel->SetSolverAllParameters(currentsolver, parameters3[0], parameters2[0], parameters1[0],
                                                       parameters22, errorMessage);

    if (work_thread.joinable())
        work_thread.join();

    progress  = 0;
    flagAbort = true;
    errorMsg.clear();
    work_thread = std::thread(&AccelModelInterface::Solve, currentProject->accelModel, currentsolver,
                              std::ref(progress), std::ref(flagAbort), std::ref(errorMsg), std::ref(status));
};

void Daizy::AddFlowAccel()
{
    currentProject->accelModel->AddFlow();
    ShowProjectTreeAccel();
};
