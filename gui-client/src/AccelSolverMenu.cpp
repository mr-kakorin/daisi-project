#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "FilesBrouses.h"
#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "GroupBoxWithTextItems.h"
#include "MiddleWidget.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>

void MiddleWidget::showAccelSolverMenu(const std::string& solver)
{
    // clearLayout(middleWidgetGrid);
    // currentsolver = solver;
    std::vector<std::string> filesIn;
    std::vector<std::string> filesOut;

    std::vector<std::string> filesInDescr;
    std::vector<std::string> filesOutDescr;

    (*currentProject)->accelModel->GetSomeSolverFileName(solver, filesIn, filesOut, filesInDescr, filesOutDescr);

    int pos = 0;
    if (filesIn.size())
    {
        std::vector<std::function<int(std::string, std::string&)>> setters(filesIn.size());
        std::vector<std::string> brouseExt(filesIn.size());
        for (int i = 0; i < filesIn.size(); i++)
        {
            brouseExt[i] = "";
            setters[i]   = std::bind(&AccelModelInterface::SetSomeSolverFileName, (*currentProject)->accelModel, solver,
                                   i, 1, std::placeholders::_1, std::placeholders::_2);
        };

        filesBrouses.push_back(new FilesBrouses("Input files"));
        filesBrouses.back()->Create(filesInDescr, brouseExt, filesIn, (*currentProject)->projectFolder, setters);
        middleWidgetGrid->addWidget(filesBrouses.back()->GetPointer(), pos, 0);
        pos++;
    }

    if (filesOut.size())
    {
        //	for (int i = 0; i < filesOut.size(); i++)
        //		filesOut[i] = GetFileName(filesOut[i]);

        groupBoxWithTextItems.push_back(new GroupBoxWithTextItems("Output files"));
        middleWidgetGrid->addWidget(groupBoxWithTextItems.back()->GetPointer(), pos, 0);
        groupBoxWithTextItems.back()->Create(filesOutDescr, filesOut);
        pos++;
    }

    std::vector<double>      params;
    std::vector<std::string> keys;

    (*currentProject)->accelModel->GetSolverParameters(solver, keys, params);

    if (params.size())
    {
        groupBoxes.push_back(new GroupBoxWithItems("Parameters"));
        middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
        groupBoxes.back()->Create(keys, params);
        pos++;
    }

    std::vector<double>                   paramsFlags;
    std::vector<std::vector<std::string>> keysFlags;

    (*currentProject)->accelModel->GetSolverParametersFlags(solver, keysFlags, paramsFlags);

    for (int i = 0; i < keysFlags.size(); i++)
    {
        groupBoxes.push_back(new GroupBoxWithItems("Parameters flags"));
        middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
        groupBoxes.back()->Create({}, {}, keysFlags[i], paramsFlags[i]);
        pos++;
    }

    QPushButton* mainMiddleWidgetButton = new QPushButton("Solve");
    mainMiddleWidgetButton->setMaximumSize(QSize(60, 40));
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AccelSolve()));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, pos, 0);
    pos++;

    AddLastFictiveBox(pos);
};