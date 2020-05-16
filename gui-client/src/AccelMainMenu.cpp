#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "FilesBrouses.h"
#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "VisualizationFunctionsAccel.h"
#include "vtkComponent.h"

#include <QtWidgets/QGridLayout>

int currentAccelSection;

void MiddleWidget::showAccelParameters(int secN)
{
    currentAccelSection = secN;

    std::vector<std::string> brouseExt;
    std::vector<std::string> brouseNames;

    (*currentProject)->accelModel->GetBrouseFlags(brouseExt, brouseNames);
    std::vector<std::string> filesNames = (*currentProject)->accelModel->GetSomeFileName();

    int pos = 0;
    if (brouseExt.size())
    {
        FilesBrouses* filesBrouses = new FilesBrouses("Parameters files");
        filesBrouses->Create(brouseNames, brouseExt, filesNames, (*currentProject)->projectFolder,
                             {std::bind(&AccelModelInterface::SetSomeParametersFromFile, (*currentProject)->accelModel,
                                        0, std::placeholders::_1, std::placeholders::_2),
                              std::bind(&AccelModelInterface::SetSomeParametersFromFile, (*currentProject)->accelModel,
                                        1, std::placeholders::_1, std::placeholders::_2)});
        middleWidgetGrid->addWidget(filesBrouses->GetPointer(), 0, 0);
        pos++;
    }

    std::vector<double>      params;
    std::vector<std::string> keys;

    std::vector<double>      paramsF;
    std::vector<std::string> keysF;

    (*currentProject)->accelModel->GetMainAccelParameters(keys, params);
    (*currentProject)->accelModel->GetMainAccelParametersFlags(keysF, paramsF);

    groupBoxes.push_back(new GroupBoxWithItems("Main parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    groupBoxes.back()->Create(keys, params, {}, 0, keysF, params);
    pos++;

    std::vector<std::vector<std::string>> keys00;
    std::vector<std::vector<double>>      p00;
    int                                   i = 0;
    (*currentProject)->accelModel->GetSectionParameters(keys00, p00);

    if (keys00.size() != 0)
    {

        QListWidget* List = new QListWidget();
        middleWidgetGrid->addWidget(List, pos, 0);
        pos++;

        List->setMaximumHeight(22 * keys00.size());

        for (int i = 0; i < keys00.size(); i++)
        {
            MyQListWidgetItem* sec = new MyQListWidgetItem();
            sec->flag              = i;
            sec->flag1             = pos;
            sec->setText("Section " + QString::number(i));
            List->addItem(sec);
        };
        connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowAccelSection(QListWidgetItem*)));

        ShowAccelSectionParams(i, pos);

        pos++;
    };

    std::vector<std::string> namesOfSequences = (*currentProject)->accelModel->GetnamesOfSequences();
    if (namesOfSequences.size())
    {
        QListWidget* ListV = new QListWidget();
        ListV->setMaximumHeight(22 * namesOfSequences.size());
        ListV->setEnabled(true);
        middleWidgetGrid->addWidget(ListV, pos, 0);
        pos++;

        for (int i = 0; i < namesOfSequences.size(); i++)
        {
            MyQListWidgetItem* secV = new MyQListWidgetItem();
            secV->flag              = i;
            secV->setText(namesOfSequences[i].c_str());
            ListV->addItem(secV);
        };
        connect(ListV, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listLinacControls(QListWidgetItem*)));
        (*VTKArray)[0]->setDataFunction(ShowLinacSequence, (*currentProject)->accelModel, 0);
        (*VTKArray)[0]->refresh(0);
    };

    if ((*currentProject)->accelModel->GetMainAccelParameterCalculated().size())
    {
        GroupBoxWithItems* MainparametersCalc = new GroupBoxWithItems("Calculated parameters");
        middleWidgetGrid->addWidget(MainparametersCalc->GetPointer(), pos, 0);
        MainparametersCalc->Create((*currentProject)->accelModel->GetMainAccelParameterNamesCalculated(),
                                   (*currentProject)->accelModel->GetMainAccelParameterCalculated(), "labels");
        pos++;
    }

    AddLastFictiveBox(pos);

    // connect(pushBrouse1, SIGNAL(clicked()), this, SLOT(fileBoundaryBrouseElements()));
};

/*void MiddleWidget::SetAccelMainParameters()
{
        bool ok = true;

        std::vector<double> params;
        for (int i = 0; i < MiddleWidgetlabelsVector.size(); i++)
        {
                params.push_back(MiddleWidgetTextEditVector[i]->toDoubleMy(&ok));
        };
        (*currentProject)->accelModel->SetMainAccelParameters(params);
        if (!ok)
        {
                QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
        }
};*/

void MiddleWidget::listShowAccelSection(QListWidgetItem* item)
{

    MyQListWidgetItem* item1 = dynamic_cast<MyQListWidgetItem*>(item);

    middleWidgetGrid->removeWidget(groupBoxes[1]->GetPointer());
    delete groupBoxes[1];
    groupBoxes.pop_back();
    currentAccelSection = item1->flag;
    ShowAccelSectionParams(item1->flag, item1->flag1);

    // showAccelParameters(item1->flag);
}

void MiddleWidget::ShowAccelSectionParams(int i, int pos)
{
    std::vector<std::vector<std::string>> keys00;
    std::vector<std::vector<double>>      p00;
    (*currentProject)->accelModel->GetSectionParameters(keys00, p00);

    groupBoxes.push_back(new GroupBoxWithItems("Accel Section " + QString::number(i)));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    groupBoxes.back()->Create(keys00[i], p00[i]);
    connect(groupBoxes.back(), SIGNAL(textChanged()), this, SLOT(SetAccelSectionParameters()));
};

void MiddleWidget::SetAccelSectionParameters()
{
    bool ok = true;

    std::vector<double> params;

    if (groupBoxes.back()->GetParameters(params))
        (*currentProject)->accelModel->SetAccelSectionParameters(currentAccelSection, params);

    (*VTKArray)[0]->refresh(0);
};

void MiddleWidget::listLinacControls(QListWidgetItem* item)
{
    MyQListWidgetItem* item1 = dynamic_cast<MyQListWidgetItem*>(item);

    (*VTKArray)[0]->setDataFunction(ShowLinacSequence, (*currentProject)->accelModel, item1->flag);
    (*VTKArray)[0]->refresh(0);
    return;
};

void MiddleWidget::showAccelParametersCalc(int secN)
{
    std::vector<std::string> Names = (*currentProject)->accelModel->GetcalcParametersNames();

    QListWidget* ListV = new QListWidget();
    ListV->setMaximumHeight(22 * Names.size());
    ListV->setEnabled(true);
    middleWidgetGrid->addWidget(ListV, 0, 0);

    for (int i = 0; i < Names.size(); i++)
    {
        MyQListWidgetItem* secV = new MyQListWidgetItem();
        secV->flag              = i;
        secV->setText(Names[i].c_str());
        ListV->addItem(secV);
    };
    connect(ListV, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(ShowAccelCalcParametersEvent(QListWidgetItem*)));

    AddLastFictiveBox(1);
};

void MiddleWidget::ShowAccelCalcParametersEvent(QListWidgetItem* item)
{
    MyQListWidgetItem* item1 = dynamic_cast<MyQListWidgetItem*>(item);

    (*VTKArray)[0]->setDataFunction(ShowAccelCalcParameters, (*currentProject)->accelModel, item1->flag);
    (*VTKArray)[0]->refresh(0);
    return;
};