#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GroupBoxWithItems.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "flagStringsResults.h"
#include "vtkComponent.h"

void MiddleWidget::showSimulationsResuts(int numberOfResultPlots)
{
    clear();

    InitTrees(numberOfResultPlots);

    for (int nr = 0; nr < numberOfResultPlots; nr++)
    {
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>> data =
            (*currentProject)->currentModel->GetDynamicsData();

        std::vector<MyTreeItem*> results(data.size());

        for (int nd = 0; nd < data.size(); nd++)
        {
            results[nd] = new MyTreeItem();
            results[nd]->setText(0, QString("Simulations results ") + QString::number(nd));
            TreeList[nr]->addTopLevelItem(results[nd]);
            results[nd]->flag1 = -1;
        }

        std::vector<std::shared_ptr<SimulationData>> simulationData =
            (*currentProject)->currentModel->GetSimulationData();

        for (int nd = 0; nd < data.size(); nd++)
        {
            MyTreeItem* timeplot = new MyTreeItem();
            timeplot->flag1      = -1;

            timeplot->setText(0, "Time plots");
            results[nd]->addChild(timeplot);

            for (int k1 = 0; k1 < simulationData[nd]->YDataFlow.size(); k1++)
            {

                MyTreeItem* Flow = new MyTreeItem();
                Flow->flag1      = -1;

                Flow->setText(0, QString("Flow ") + QString::number(k1));
                timeplot->addChild(Flow);

                for (int k2 = 0; k2 < simulationData[nd]->YDataFlow[k1].size(); k2++)
                {

                    itemvector[nr].push_back(new MyTreeItem());
                    itemvector[nr].back()->flag  = "Time";
                    itemvector[nr].back()->flag1 = k1;
                    itemvector[nr].back()->flag2 = k2;
                    itemvector[nr].back()->flag4 = 0;
                    itemvector[nr].back()->flag5 = nr;
                    itemvector[nr].back()->flag6 = nd;
                    itemvector[nr].back()->flag3 = simulationData[nd]->dataFlags[0][k2];

                    itemvector[nr].back()->setText(0, simulationData[nd]->dataFlags[0][k2].c_str());
                    Flow->addChild(itemvector[nr].back());
                }
            }

            for (int k1 = 0; k1 < simulationData[nd]->YDataElectrode.size(); k1++)
            {

                MyTreeItem* Flow = new MyTreeItem();
                Flow->flag1      = -1;

                Flow->setText(0, QString("Electrode ") + QString::number(k1));
                timeplot->addChild(Flow);

                for (int k2 = 0; k2 < simulationData[nd]->YDataElectrode[k1].size(); k2++)
                {

                    itemvector[nr].push_back(new MyTreeItem());
                    itemvector[nr].back()->flag  = "Time";
                    itemvector[nr].back()->flag1 = k1;
                    itemvector[nr].back()->flag2 = k2;
                    itemvector[nr].back()->flag4 = 1;
                    itemvector[nr].back()->flag5 = nr;
                    itemvector[nr].back()->flag6 = nd;
                    itemvector[nr].back()->flag3 = simulationData[nd]->dataFlags[1][k2];

                    itemvector[nr].back()->setText(0, simulationData[nd]->dataFlags[1][k2].c_str());
                    Flow->addChild(itemvector[nr].back());
                }
            }
        }

        std::vector<std::string> currentNames;

        for (int t = 0; t < data.size(); t++)
        {
            for (int t1 = 0; t1 < data[t].size(); t1++)
            {
                if (data[t][t1][0]->StartTime() == -1)
                    continue;

                MyTreeItem* time = new MyTreeItem();
                time->flag1      = -1;

                time->setText(0,
                              (std::string("start time, ns ") + std::to_string(data[t][t1][0]->StartTime())).c_str());
                results[t]->addChild(time);
                std::vector<int> ParticlesFlowsTypes = (*currentProject)->currentModel->GetNumberParticlesFlowsTypes();
                for (int k = 0; k < (*currentProject)->currentModel->GetNumberParticlesFlows(); k++)
                {
                    if (ParticlesFlowsTypes[k] < 5)
                    {
                        switch ((*currentProject)->problemType)
                        {
                        case 1:
                            currentNames = flagStrings::ResultsNames2d;
                            break;
                        case 2:
                            currentNames = flagStrings::ResultsNames2daxs;
                            break;
                        case 3:
                            currentNames = flagStrings::ResultsNames2dpolar;
                            break;
                        case 4:
                            currentNames = flagStrings::ResultsNames3d;
                            break;
                        }
                    }
                    else
                    {
                        switch ((*currentProject)->problemType)
                        {
                        case 1:
                            currentNames = flagStrings::ResultsNames2dLinac;
                            break;
                        case 2:
                            currentNames = flagStrings::ResultsNames2daxsLinac;
                            break;
                        case 3:
                            currentNames = flagStrings::ResultsNames2dpolar;
                            break;
                        case 4:
                            currentNames = flagStrings::ResultsNames3d;
                            break;
                        }
                    };

                    MyTreeItem* flow = new MyTreeItem();
                    flow->flag1      = -1;

                    flow->setText(0, (std::string("flow") + std::to_string(k)).c_str());
                    time->addChild(flow);

                    for (int i = 0; i < currentNames.size(); i++)
                    {
                        itemvector[nr].push_back(new MyTreeItem());
                        itemvector[nr].back()->flag3 = "trace";
                        itemvector[nr].back()->flag  = currentNames[i];
                        itemvector[nr].back()->setText(0, itemvector[nr].back()->flag.c_str());
                        itemvector[nr].back()->flag1 = t;
                        itemvector[nr].back()->flag2 = t1;
                        itemvector[nr].back()->flag4 = k;
                        itemvector[nr].back()->flag5 = nr;

                        flow->addChild(itemvector[nr].back());
                    };
                };
            }
        }
        connect(TreeList[nr], SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
                SLOT(listShowResultClick(QTreeWidgetItem*, int)));
    }

    groupBoxes.push_back(new GroupBoxWithItems("Main parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), numberOfResultPlots, 0);
    groupBoxes.back()->Create({"Number of particles", "L end (-1 means Lmax)"}, std::vector<int>{100, -1});
}

void MiddleWidget::listShowResultClick(QTreeWidgetItem* item, int)
{
    std::vector<std::shared_ptr<SimulationData>> simulationData = (*currentProject)->currentModel->GetSimulationData();

    MyTreeItem* item1 = dynamic_cast<MyTreeItem*>(item);

    if (-1 == item1->flag1)
        return;

    std::vector<int>    props = {-1, -1, -1};
    std::vector<double> props1;
    bool                ok = groupBoxes[0]->GetParameters(props1);

    int numberOfTraces = props1[0];

    int Distrtype = 0;

    if (item1->flag == std::string("Time"))
    {
        (*VTKArray)[item1->flag5]->setDataFunction(ShowSimulationDataPlot, (*currentProject)->currentModel,
                                                   item1->flag6, item1->flag4, item1->flag1, item1->flag2,
                                                   item1->flag3);
        (*VTKArray)[item1->flag5]->refresh(0);
        return;
    };

    std::shared_ptr<DynamicsData> data =
        (*currentProject)->currentModel->GetDynamicsData()[item1->flag1][item1->flag2][item1->flag4];

    switch ((*currentProject)->problemType)
    {
    case 1:
        if (item1->flag == "X(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("X, m"), 0,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Y(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("Y, m"), 1,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Z(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("Z, m"), 2,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PX(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("PX(t)"), 3,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PY(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("PY(t)"), 4,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PX(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("PX(t)"), 3,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PY(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("PY(t)"), 4,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "X(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("X, m"), 0,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Y(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("Y, m"), 1,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Charge(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("Charge, cl"),
                6, item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Energy(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy, (*currentProject)->currentModel, data,
                                                       std::string("time, ns"), item1->flag4, numberOfTraces, Distrtype,
                                                       props, 0);

        if (item1->flag == "Energy(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy, (*currentProject)->currentModel, data,
                                                       std::string("Z, m"), item1->flag4, numberOfTraces, Distrtype,
                                                       props, 0);

        if (item1->flag == "Traces in XY plane")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane, (*currentProject)->currentModel,
                                                       data, item1->flag4, std::string("X, m"), std::string("Y, m"),
                                                       (*currentProject)->problemType, numberOfTraces, 1, Distrtype,
                                                       props);

        break;
    case 2:

        if (item1->flag == "R(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("R, m"), 0,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "R(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("R, m"), 0,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Z(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("Z, m"), 1,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Traces in RZ plane")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane, (*currentProject)->currentModel,
                                                       data, item1->flag4, std::string("R, m"), std::string("Z, m"),
                                                       (*currentProject)->problemType, numberOfTraces, 0, Distrtype,
                                                       props);

        if (item1->flag == "Traces in ZR plane")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane, (*currentProject)->currentModel,
                                                       data, item1->flag4, std::string("Z, m"), std::string("R, m"),
                                                       (*currentProject)->problemType, numberOfTraces, 1, Distrtype,
                                                       props);

        if (item1->flag == "PR(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("PR(t)"), 3,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PR(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("PR(t)"), 3,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PZ(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("Z, m"), std::string("PZ(t)"), 4,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "PZ(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("PZ(t)"), 4,
                item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Charge(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowTracesT, (*currentProject)->currentModel, data, std::string("time, ns"), std::string("Charge(t)"),
                6, item1->flag4, numberOfTraces, Distrtype, props, (*currentProject)->problemType);

        if (item1->flag == "Energy(t)")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy, (*currentProject)->currentModel, data,
                                                       std::string("time, ns"), item1->flag4, numberOfTraces, Distrtype,
                                                       props, 0);

        if (item1->flag == "Energy(Z)")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy, (*currentProject)->currentModel, data,
                                                       std::string("Z, m"), item1->flag4, numberOfTraces, Distrtype,
                                                       props, 0);

        if (item1->flag == "Traces in XY plane")
            (*VTKArray)[item1->flag5]->setDataFunction(ShowXY, (*currentProject)->currentModel, data, item1->flag4,
                                                       numberOfTraces, Distrtype, props);

        /*if (item1->flag == flagStrings::ResultsNames2daxs[0])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane,
        (*currentProject)->currentModel, data, item1->flag4, std::string("R, m"), std::string("Z, m"),
        (*currentProject)->problemType, numberOfTraces, 0, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        if (item1->flag == flagStrings::ResultsNames2daxs[1])
        {
        //	(*VTKArray)[item1->flag5]->setDataFunction(ShowEnergyEZ, (*currentProject)->currentModel, data,
        item1->flag4, numberOfTraces, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        if (item1->flag == flagStrings::ResultsNames2daxs[2])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy, (*currentProject)->currentModel, data,
        item1->flag4, numberOfTraces, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        if (item1->flag == flagStrings::ResultsNames2daxs[3])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergyEphiPolar, (*currentProject)->currentModel, data,
        item1->flag4, numberOfTraces, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };


        if (item1->flag == flagStrings::ResultsNames2daxs[4])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(ShowXY, (*currentProject)->currentModel, data, item1->flag4,
        numberOfTraces, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        if (item1->flag == flagStrings::ResultsNames2daxs[5])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane,
        (*currentProject)->currentModel, data, item1->flag4, std::string("R, m"), std::string("Z, m"),
        (*currentProject)->problemType, numberOfTraces, 1, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        if (item1->flag == flagStrings::ResultsNames2daxs[6])
        {
                (*VTKArray)[item1->flag5]->setDataFunction(Show3PositionsWithGeometryPlane,
        (*currentProject)->currentModel, data, item1->flag4, std::string("R, m"), std::string("Z, m"),
        (*currentProject)->problemType, numberOfTraces, 1, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
        };

        /*for (int i = 7; i < flagStrings::ResultsNames2daxs.size(); i++)
        {
                if (item1->flag == flagStrings::ResultsNames2daxs[i])
                {
                        (*VTKArray)[item1->flag5]->setDataFunction(Show2daxsTrace, (*currentProject)->currentModel,
        data, item1->flag4, numberOfTraces, i - 7, Distrtype, props);
                        (*VTKArray)[item1->flag5]->refresh(0);
                        return;
                };
        }*/
        break;

    case 3:
        if (item1->flag == flagStrings::ResultsNames2dpolar[0])
        {
            (*VTKArray)[item1->flag5]->setDataFunction(ShowPositionsWithGeometryPlane, (*currentProject)->currentModel,
                                                       data, item1->flag4, std::string("X, m"), std::string("Y, m"),
                                                       (*currentProject)->problemType, numberOfTraces, 0, Distrtype,
                                                       props);
            (*VTKArray)[item1->flag5]->refresh(0);
            return;
        };
        if (item1->flag == flagStrings::ResultsNames2dpolar[4])
        {
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergyPolar, (*currentProject)->currentModel, data,
                                                       item1->flag4, numberOfTraces, Distrtype, props);
            (*VTKArray)[item1->flag5]->refresh(0);
            return;
        };
        if (item1->flag == flagStrings::ResultsNames2dpolar[5])
        {
            (*VTKArray)[item1->flag5]->setDataFunction(
                ShowPositionsWithGeometryPlaneSpecial, (*currentProject)->currentModel, data, item1->flag4,
                std::string("X, m"), std::string("Y, m"), (*currentProject)->problemType, numberOfTraces, Distrtype,
                props);
            (*VTKArray)[item1->flag5]->refresh(0);
            return;
        };

        if (item1->flag == flagStrings::ResultsNames2dpolar[6])
        {
            (*VTKArray)[item1->flag5]->setDataFunction(ShowEnergyEphiPolar, (*currentProject)->currentModel, data,
                                                       item1->flag4, numberOfTraces, Distrtype, props);
            (*VTKArray)[item1->flag5]->refresh(0);
            return;
        };

        break;

    case 4:
        for (int i = 0; i < 6; i++)
        {
            if (item1->flag == flagStrings::ResultsNames3d[i])
            {
                (*VTKArray)[item1->flag5]->setDataFunction(Show3dTrace, (*currentProject)->currentModel, data,
                                                           item1->flag4, numberOfTraces, i, Distrtype, props);
                (*VTKArray)[item1->flag5]->refresh(0);
                return;
            };
        }
        if (item1->flag == flagStrings::ResultsNames3d[6])
        {
            //	(*VTKArray)[item1->flag5]->setDataFunction(ShowEnergy3d, (*currentProject)->currentModel, data,
            // item1->flag4, numberOfTraces, Distrtype, props);
            (*VTKArray)[item1->flag5]->refresh(0);
            return;
        };
        break;
    }
    (*VTKArray)[item1->flag5]->refresh(0);
};