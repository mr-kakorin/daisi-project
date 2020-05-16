#include "FilesBrouses.h"
#include "FlagStringsD.h"
#include "GroupBoxWithItems.h"
#include "MiddleWidget.h"

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "MyTreeItem.h"
#include "flagStringsTree.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>

void MiddleWidget::ShowAddBoundary()
{
    std::vector<std::function<int(std::string, std::string&)>> setters(2);
    std::vector<std::string> brouseNames = {"Input *.dat file", "Input config file"};
    std::vector<std::string> brouseExt   = {"*.dat", ""};
    std::vector<std::string> filesIn     = {"", ""};
    setters[0] = std::bind(&MiddleWidget::AddBoundaryEvent, this, std::placeholders::_1, std::placeholders::_2);
    setters[1] = std::bind(&MiddleWidget::AddBoundariesEvent, this, std::placeholders::_1, std::placeholders::_2);

    filesBrouses.push_back(new FilesBrouses("Input files"));
    filesBrouses.back()->Create(brouseNames, brouseExt, filesIn, (*currentProject)->projectFolder, setters);
    middleWidgetGrid->addWidget(filesBrouses.back()->GetPointer(), 0, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Change geometry"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);
    groupBoxes.back()->Create({"Changed point", "Changing point"}, std::vector<double>{0, 0});

    mainMiddleWidgetButton = new QPushButton(tr("Change geom"));
    mainMiddleWidgetButton->setMaximumSize(QSize(150, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 2, 0);
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ChangeGeom()));

    AddLastFictiveBox(3);
};
int MiddleWidget::AddBoundaryEvent(const std::string& input, std::string& errorMsg)
{
    (*currentProject)->currentModel->AddBoundary(input, errorMsg);

    MyTreeItem* newItem = new MyTreeItem();
    int         i       = (*currentProject)->currentModel->GetNumberBoundaries();
    newItem->setText(
        0, QApplication::translate("MainWindow", (std::string("boundary") + std::to_string(i - 1)).c_str(), 0));
    newItem->flag  = flagStrings::boundariesList;
    newItem->flag1 = i - 1;
    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);
    newItem->setSelected(true);

    emit ShowProjectTreeSim();
    return 1;
};

int MiddleWidget::AddBoundariesEvent(const std::string& input, std::string& errorMsg)
{
    (*currentProject)->currentModel->AddBoundaries(input, errorMsg);

    emit ShowProjectTreeSim();
    return 1;
};

void MiddleWidget::ChangeGeom()
{
    std::vector<double> parameters;
    if (groupBoxes[0]->GetParameters(parameters))
        (*currentProject)->currentModel->ChangeGeom(parameters[0], parameters[1]);
    else
        QMessageBox::critical(this, "Daisi warning", "Incorrect input");
};