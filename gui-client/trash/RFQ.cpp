#include "DaisiClient.h"
#include "FlagStringsD.h"
#include "ProjectTree.h"
#include "VizualizationFunctions.h"
#include "vtkComponent.h"

void Daizy::showRFQMainParameters()
{
    /*	clearLayout(middleWidgetGrid);

            std::vector<double> params = currentProject->currentModel->GetLinacParameters();
    //	std::vector<std::vector<double>> paramsAppr = currentProject->currentModel->GetLinacApproximationParameters();

            QGridLayout *vbox1 = new QGridLayout;
            QGroupBox *groupBox1 = new QGroupBox(tr("RFQ main parameters"));

            middleWidgetGrid->addWidget((QWidget*)groupBox1, 0, 0);

            groupBox1->setLayout(vbox1);
            groupBox1->setMaximumHeight(500);

            MiddleWidgetTextEditVector.clear();
            for (int j = 0; j < flagStrings::RFQMainParameters.size(); j++)
            {
                    QLabel* label = new QLabel();
                    label->setText(flagStrings::RFQMainParameters[j].c_str());
                    vbox1->addWidget(label, j, 0);

                    MiddleWidgetTextEditVector.push_back(new MyQTextEdit());
                    MiddleWidgetTextEditVector.back()->setMaximumSize(QSize(16777215, 25));
                    MiddleWidgetTextEditVector.back()->setText(QString::number(params[j]));
                    vbox1->addWidget(MiddleWidgetTextEditVector.back(), j, 1);
            };

            std::vector<QGridLayout *> vboxSections(paramsAppr.size());
            std::vector<QGroupBox *> groupBoxSections(paramsAppr.size());
            MiddleWidgetTextEditDoubleVector.resize(paramsAppr.size());

            for (int section = 0; section < paramsAppr.size(); section++)
            {
                    vboxSections[section] = new QGridLayout;
                    groupBoxSections[section] = new QGroupBox(flagStrings::RFQSectionsNames[section].c_str());
                    MiddleWidgetTextEditDoubleVector[section].resize(paramsAppr[section].size());
                    groupBoxSections[section]->setLayout(vboxSections[section]);

                    middleWidgetGrid->addWidget((QWidget*)groupBoxSections[section], section + 1, 0);

                    for (int j = 0; j < paramsAppr[section].size(); j++)
                    {
                            QLabel* label = new QLabel();
                            label->setText(flagStrings::RFQSectionsNamesFields[j].c_str());

                            if (section > paramsAppr.size() - 3)
                                    label->setText(flagStrings::RFQMatcherNames[j].c_str());

                            int pos = 0;
                            int pos0 = j;

                            if (j == 2 || j == 4)
                                    pos = 2;

                            if (j == 2 || j == 1)
                                    pos0 = 1;

                            if (j == 4 || j == 3)
                                    pos0 = 2;

                            vboxSections[section]->addWidget(label, pos0, pos);

                            MiddleWidgetTextEditDoubleVector[section][j] = new MyQTextEdit();
                            MiddleWidgetTextEditDoubleVector[section][j]->setMaximumSize(QSize(16777215, 25));
                            MiddleWidgetTextEditDoubleVector[section][j]->setText(QString::number(paramsAppr[section][j]));
                            vboxSections[section]->addWidget(MiddleWidgetTextEditDoubleVector[section][j], pos0, pos +
    1);
                    };
            };


            MiddleWidgetButton1 = new QPushButton(tr("Set main parameters"));
            MiddleWidgetButton1->setMaximumSize(QSize(150, 30));
            middleWidgetGrid->addWidget(MiddleWidgetButton1, paramsAppr.size()+1,0);
            connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(SetLinacParameters()));



            QListWidget* List = new QListWidget();

            for (int i = 0; i < flagStrings::RFQControlsNames.size(); i++)
            {
                    QListWidgetItem* newitem = new QListWidgetItem();
                    newitem->setText(flagStrings::RFQControlsNames[i].c_str());
                    List->addItem(newitem);
            };

            middleWidgetGrid->addWidget(List, paramsAppr.size() + 2, 0, 0);

            connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listLinacControls(QListWidgetItem*)));*/
}

void Daizy::SetLinacParameters()
{
    bool ok = true;

    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }

    std::vector<std::vector<double>> paramsAppr(MiddleWidgetTextEditDoubleVector.size());
    for (int j = 0; j < MiddleWidgetTextEditDoubleVector.size(); j++)
    {
        paramsAppr[j].resize(MiddleWidgetTextEditDoubleVector[j].size());
        for (int i = 0; i < MiddleWidgetTextEditDoubleVector[j].size(); i++)
        {
            paramsAppr[j][i] = MiddleWidgetTextEditDoubleVector[j][i]->toDoubleMy(&ok);
        }
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    //	currentProject->currentModel->SetLinacParameters(params, paramsAppr, currentProject->projectFolder);
}
void Daizy::SetSyncPhases()
{
    //	currentProject->currentModel->SetControlFunction(MiddleWidgetTextEdit1->toPlainText().toStdString(), 1);
}
void Daizy::SetAccelEff()
{
    //	currentProject->currentModel->SetControlFunction(MiddleWidgetTextEdit1->toPlainText().toStdString(), 0);
}

/*void Daizy::ShowFlowAccel(int flow)
{
        clearLayout(middleWidgetGrid);

        std::vector<double> params = currentProject->currentModel->GetLinacFlowParameters(flow);

        
        QGridLayout *vbox1 = new QGridLayout;
        QGroupBox *groupBox1 = new QGroupBox(tr("RFQ main parameters"));

        middleWidgetGrid->addWidget((QWidget*)groupBox1, 0, 0);

        groupBox1->setLayout(vbox1);
        groupBox1->setMaximumHeight(600);

        MiddleWidgetTextEditVector.clear();

        int s = 0;

        for (int j = 0; j < flagStrings::RFQFlowParameters.size(); j++)
        {
                QLabel* label = new QLabel();
                label->setText(flagStrings::RFQFlowParameters[j].c_str());
                vbox1->addWidget(label, j, 0);

                MiddleWidgetTextEditVector.push_back(new MyQTextEdit());
                MiddleWidgetTextEditVector.back()->setMaximumSize(QSize(16777215, 25));
                MiddleWidgetTextEditVector.back()->setText(QString::number(params[j]));
                vbox1->addWidget(MiddleWidgetTextEditVector.back(), j, 1);
                s++;
        };

        MiddleWidgetButton1 = new QPushButton(tr("Set flow parameters"));
        MiddleWidgetButton1->setMaximumSize(QSize(150, 30));
        vbox1->addWidget(MiddleWidgetButton1, s, 0);
        currentFlow = flow;
        connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(SetLinacFlowParameters()));
        s++;

        MiddleWidgetButton2 = new QPushButton(tr("Generate RFQ for this flow"));
        MiddleWidgetButton2->setMaximumSize(QSize(150, 30));
        vbox1->addWidget(MiddleWidgetButton2, s, 0);
        currentFlow = flow;
        connect(MiddleWidgetButton2, SIGNAL(clicked()), this, SLOT(GenerateRFQForFlow()));
        s++;


        MiddleWidgetButton3 = new QPushButton(tr("Simulate dynamics"));
        MiddleWidgetButton3->setMaximumSize(QSize(150, 30));
        vbox1->addWidget(MiddleWidgetButton3, s, 0);
        currentFlow = flow;
        connect(MiddleWidgetButton3, SIGNAL(clicked()), this, SLOT(SimulateDynRFQ()));
        s++;


        MiddleWidgetButton4 = new QPushButton(tr("Calculate acceptance"));
        MiddleWidgetButton4->setMaximumSize(QSize(150, 30));
        vbox1->addWidget(MiddleWidgetButton4, s, 0);
        currentFlow = flow;
        connect(MiddleWidgetButton4, SIGNAL(clicked()), this, SLOT(CalculateRFQAcc()));
        s++;


        QLabel* label = new QLabel();
        label->setText("file");
        vbox1->addWidget(label, s, 0);

        MiddleWidgetTextEdit1 = new MyQTextEdit();
        MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));
        vbox1->addWidget(MiddleWidgetTextEdit1, s, 1);
        s++;


        

        MiddleWidgetButton5 = new QPushButton(tr("Export to LIDOS"));
        MiddleWidgetButton5->setMaximumSize(QSize(150, 30));
        vbox1->addWidget(MiddleWidgetButton5, s, 0);
        connect(MiddleWidgetButton5, SIGNAL(clicked()), this, SLOT(ExportToLidos()));
        s++;

        currentFlow = flow;
}*/
void Daizy::CalculateRFQAcc()
{
    if (work_thread.joinable())
        work_thread.join();

    //	work_thread = std::thread(&ModelInterface::CalculateRFQAcc, currentProject->currentModel, currentFlow);
};
void Daizy::ExportToLidos()
{
    std::string filename =
        currentProject->projectFolder + "\\" + MiddleWidgetTextEdit1->toPlainText().toStdString() + ".dat";
    //	currentProject->currentModel->ExportToLidos(filename, currentFlow);
};
void Daizy::SimulateDynRFQ(){
    //	currentProject->currentModel->SimulateLinacFlowDynamics(currentFlow);
};
void Daizy::SetLinacFlowParameters()
{
    bool ok = true;

    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    //	currentProject->currentModel->SetLinacFlowParameters(currentFlow, params);
};

void Daizy::GenerateRFQForFlow()
{
    int succes;
    //	currentProject->currentModel->GenerateRFQForFlow(succes, currentFlow);

    if (succes == 0)
    {
        QMessageBox messageBox;
        messageBox.critical(0, "Error", "It is impossible to generate RFQ for this channel aperture");
        messageBox.setFixedSize(500, 200);

        messageBox.exec();
    }
};

void Daizy::ShowRFQCavityParameters()
{
    clearLayout(middleWidgetGrid);

    QListWidget* List = new QListWidget();

    for (int i = 0; i < flagStrings::RFQCavityParametersNames.size(); i++)
    {
        QListWidgetItem* newitem = new QListWidgetItem();
        newitem->setText(flagStrings::RFQCavityParametersNames[i].c_str());
        List->addItem(newitem);
    };

    middleWidgetGrid->addWidget(List, 2, 0);

    connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowRFQCavityParameters(QListWidgetItem*)));
}

void Daizy::listShowRFQCavityParameters(QListWidgetItem* item)
{
    std::string name = item->text().toStdString();

    for (int i = 0; i < flagStrings::RFQCavityParametersNames.size(); i++)
    {
        if (name == flagStrings::RFQCavityParametersNames[i])
        {
            VTKArray[0]->setDataFunction(ShowRFQRFQCavityParametersPlots, currentProject->currentModel, i);
            VTKArray[0]->refresh(0);
            return;
        };
    }
};

void Daizy::ShowRFQSimParams()
{
    clearLayout(middleWidgetGrid);

    std::vector<double> params = currentProject->currentModel->GetLinacDynSimParameters();

    QGridLayout* vbox1     = new QGridLayout;
    QGroupBox*   groupBox1 = new QGroupBox(tr("Dyn sim parameters"));

    middleWidgetGrid->addWidget((QWidget*)groupBox1, 0, 0);

    groupBox1->setLayout(vbox1);
    groupBox1->setMaximumHeight(300);

    MiddleWidgetTextEditVector.clear();
    for (int j = 0; j < flagStrings::LinacSimParameters.size(); j++)
    {
        QLabel* label = new QLabel();
        label->setText(flagStrings::LinacSimParameters[j].c_str());
        vbox1->addWidget(label, j, 0);

        MiddleWidgetTextEditVector.push_back(new MyQTextEdit());
        MiddleWidgetTextEditVector.back()->setMaximumSize(QSize(16777215, 25));
        MiddleWidgetTextEditVector.back()->setText(QString::number(params[j]));
        vbox1->addWidget(MiddleWidgetTextEditVector.back(), j, 1);
    };

    MiddleWidgetButton1 = new QPushButton(tr("Set parameters"));
    MiddleWidgetButton1->setMaximumSize(QSize(150, 30));
    vbox1->addWidget(MiddleWidgetButton1);
    connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(SetLinacDynSimParameters()));
}

void Daizy::SetLinacDynSimParameters()
{
    bool                ok = true;
    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    //	currentProject->currentModel->SetLinacDynSimParameters(params);
};

void Daizy::RFQOpt()
{
    bool ok = true;

    clearLayout(middleWidgetGrid);

    QGridLayout* vbox1     = new QGridLayout;
    QGroupBox*   groupBox1 = new QGroupBox(tr("RFQ main parameters"));

    middleWidgetGrid->addWidget((QWidget*)groupBox1, 0, 0);

    groupBox1->setLayout(vbox1);
    groupBox1->setMaximumHeight(100);

    std::vector<std::string> optParams = {"Number of threads", "Number of Agents"};

    std::vector<double> params = {4.0, 200.0};

    MiddleWidgetTextEditVector.clear();
    for (int j = 0; j < optParams.size(); j++)
    {
        QLabel* label = new QLabel();
        label->setText(optParams[j].c_str());
        vbox1->addWidget(label, j, 0);

        MiddleWidgetTextEditVector.push_back(new MyQTextEdit());
        MiddleWidgetTextEditVector.back()->setMaximumSize(QSize(16777215, 25));
        MiddleWidgetTextEditVector.back()->setText(QString::number(params[j]));
        vbox1->addWidget(MiddleWidgetTextEditVector.back(), j, 1);
    };

    middleWidgetGrid->addWidget(groupBox1, 0, 0);

    mainMiddleWidgetButton = new QPushButton(tr("RFQ initial Opt"));
    mainMiddleWidgetButton->setMaximumSize(QSize(170, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 1, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(RFQinitialOpt()));

    MiddleWidgetButton1 = new QPushButton(tr("RFQ output Em Opt"));
    MiddleWidgetButton1->setMaximumSize(QSize(170, 30));
    middleWidgetGrid->addWidget(MiddleWidgetButton1, 2, 0);

    connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(RFQEmOpt()));

    MiddleWidgetButton2 = new QPushButton(tr("RFQ max acc"));
    MiddleWidgetButton2->setMaximumSize(QSize(170, 30));
    middleWidgetGrid->addWidget(MiddleWidgetButton2, 3, 0);

    connect(MiddleWidgetButton2, SIGNAL(clicked()), this, SLOT(RFQAccOpt()));

    MiddleWidgetButton3 = new QPushButton(tr("RFQ max acc"));
    MiddleWidgetButton3->setMaximumSize(QSize(170, 30));
    middleWidgetGrid->addWidget(MiddleWidgetButton3, 4, 0);

    connect(MiddleWidgetButton3, SIGNAL(clicked()), this, SLOT(RFQMatcherOpt()));
};

void Daizy::RFQMatcherOpt()
{
    bool ok = true;

    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    if (work_thread.joinable())
        work_thread.join();

    flagAbort = true;
    progress  = 0;
    //	work_thread = std::thread(&ModelInterface::RFQMatcherOpt, currentProject->currentModel, std::ref(flagAbort),
    // std::ref(progress), params);
};

void Daizy::RFQAccOpt()
{
    bool ok = true;

    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    if (work_thread.joinable())
        work_thread.join();

    flagAbort = true;
    progress  = 0;
    //	work_thread = std::thread(&ModelInterface::RFQAccOpt, currentProject->currentModel, std::ref(flagAbort),
    // std::ref(progress), params);
};

void Daizy::RFQEmOpt()
{
    bool ok = true;

    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    if (work_thread.joinable())
        work_thread.join();

    flagAbort = true;
    progress  = 0;
    //	work_thread = std::thread(&ModelInterface::RFQEmOpt, currentProject->currentModel, std::ref(flagAbort),
    // std::ref(progress), params);
};
void Daizy::RFQinitialOpt()
{
    bool                ok = true;
    std::vector<double> params;
    for (int j = 0; j < MiddleWidgetTextEditVector.size(); j++)
    {
        params.push_back(MiddleWidgetTextEditVector[j]->toDoubleMy(&ok));
    }
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    if (work_thread.joinable())
        work_thread.join();

    flagAbort = true;
    progress  = 0;
    //	work_thread = std::thread(&ModelInterface::SearchStartParameters, currentProject->currentModel,
    // std::ref(flagAbort),  std::ref(progress), params);
};