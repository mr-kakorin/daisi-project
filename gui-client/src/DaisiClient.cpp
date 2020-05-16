#include <thread>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>


#include "DaisiClient.h"
#include "vtkComponent.h"
#ifdef NUCL
#include "namesNucl.h"
#else
#include "names.h"
#endif
#include "Menu.h"
#include "MiddleWidget.h"

#include "GeneralTools.h"
#include <QVTKWidget.h>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTreeWidgetItem>
#include <QtWidgets/QVBoxLayout>

QString statusLast;

Daizy::Daizy(QWidget* parent) : QMainWindow(parent)
{
    currentProject      = NULL;
    numberOfResultPlots = 1;
    flagExit            = true;
    flagAbort           = false;

    progress = 0;
    if (this->objectName().isEmpty())
        this->setObjectName("MainWindow");
    this->resize(1280, 720);
    this->setWindowModality(Qt::WindowModal);

    centralWidget = new QWidget(this);
    centralWidget->setObjectName("centralWidget");
    this->setWindowTitle(codeName.c_str());
    // this->setWindowIcon();
    menuBar = new QMenuBar(this);
    menuBar->setGeometry(QRect(0, 0, 1286, 21));
    this->setMenuBar(menuBar);

    timerViz = new QTimer(this);

    // contextMenu = new QMenu(this);
    //	contextMenu->addAction(new QAction("Action 1", this));

    splitter = new QSplitter(centralWidget);
    splitter->setObjectName(QStringLiteral("splitter"));
    splitter->setOrientation(Qt::Horizontal);

    // splitterStatus = new QSplitter(leftWidget);
    //	splitterStatus->setObjectName(QStringLiteral("splitter"));
    // splitterStatus->setOrientation(Qt::Vertical);

    leftWidget = new QWidget(splitter);
    leftWidget->setObjectName(QStringLiteral("Middle"));
    leftWidgetGrid = new QGridLayout(leftWidget);
    leftWidgetTree = new QTreeWidget();
    leftWidgetTree->setObjectName(QStringLiteral("Left"));

    QTreeWidgetItem* header = new QTreeWidgetItem();
    header->setText(0, "Model builder");
    leftWidgetTree->setHeaderItem(header);

    leftWidgetGrid->addWidget(leftWidgetTree, 0, 0);

    progressBarLocal = new QProgressBar(leftWidget);
    progressBarLocal->setObjectName(QStringLiteral("progressBar"));
    progressBarLocal->setValue(0);
    leftWidgetGrid->addWidget(progressBarLocal, 1, 0);

    progressBar = new QProgressBar(leftWidget);
    progressBar->setObjectName(QStringLiteral("progressBar"));
    progressBar->setValue(0);
    leftWidgetGrid->addWidget(progressBar, 2, 0);

    monitor = new QTextEdit(leftWidget);
    monitor->setMaximumHeight(80);
    leftWidgetGrid->addWidget(monitor, 3, 0);
    scrollBar = monitor->verticalScrollBar();

    buttonAbort = new QPushButton(tr("Abort simulations"));
    leftWidgetGrid->addWidget(buttonAbort, 5, 0);

    connect(buttonAbort, SIGNAL(clicked()), this, SLOT(AbortSimulations()));

    labelTime = new QLabel();
    labelTime->setText("Time");
    leftWidgetGrid->addWidget(labelTime, 4, 0);

    progressBarThread = std::thread(&Daizy::setProgress, this);

    splitter->addWidget(leftWidget);

    middleWidget = new QWidget(splitter);
    middleWidget->setObjectName("Middle");
    splitter->addWidget(middleWidget);
    // middleWidgetGrid = new QGridLayout(middleWidget);
    middleWidgetAggregator = new MiddleWidget(middleWidget, &currentProject, &VTKArray);

    rigthWidget = new QWidget(splitter);
    rigthWidget->setObjectName("Rigth");
    splitter->addWidget(rigthWidget);

    rightLayout = new QVBoxLayout(rigthWidget);

    rigthWidgetGrid = new QGridLayout();

    rightLayout->addLayout(rigthWidgetGrid);

    graphicsArray.push_back(new QVTKWidget());
    graphicsArray[0]->setMinimumSize(QSize(0, 0));
    rigthWidgetGrid->addWidget(graphicsArray[0], 0, 0);
    VTKArray.push_back(new vtkComponent());
    VTKArray[0]->setWidget(graphicsArray.back());

    /*SolutionParamsGrid = new QGridLayout;
    SolutionParamsBox = new QGroupBox("Solution parameters");
    SolutionParamsBox->setMaximumHeight(120);
    SolutionParamsBox->setLayout(SolutionParamsGrid);

    rightLayout->addWidget((QWidget*)SolutionParamsBox, 1, 0);*/

    splitter->setGeometry(QRect(5, 5, 1270, 690));

    QList<int> sizes = splitter->sizes();
    sizes[0]         = 250;
    sizes[1]         = 250;
    sizes[2]         = 770;
    splitter->setSizes(sizes);

    horizontalLayout = new QHBoxLayout(centralWidget);
    horizontalLayout->setSpacing(6);
    horizontalLayout->setContentsMargins(11, 11, 11, 11);
    horizontalLayout->setObjectName("horizontalLayout");
    horizontalLayout->addWidget(splitter);

    middleWidgetAggregator->ShowCreateNewProject();
    this->setCentralWidget(centralWidget);
    connect(this, SIGNAL(showMBSignal()), this, SLOT(ShowMB()));

    mainMenu = new MainMenu();
    mainMenu->MenuShow(menuBar, this, &currentProject, &VTKArray, middleWidgetAggregator->GetMiddleWidgetGrid());
    MyConnect();
    QMetaObject::connectSlotsByName(this);

    // regressionProcessor();
    // middleWidgetAggregator->ShowCreateNewProject();

    /*	std::string path = "../Projects/test/test.dproj";


            if (!(currentProject))
                    (currentProject) = new Dproject::project();

            std::string errorMsgOpen;

            (currentProject)->LoadProject(path, version, errorMsgOpen);

            if (errorMsgOpen.size())
            {
                    QMessageBox::critical(this, "Daisi error", errorMsgOpen.c_str());
                    return;
            };


            if ((currentProject)->problemType >= 7)
                    ShowProjectTreeAccel();
            else
                    ShowProjectTreeSim();

            middleWidgetAggregator->clear();
            middleWidgetAggregator->showSummary();*/
}

Daizy::~Daizy()
{
    if (work_thread.joinable())
        work_thread.join();

    if (Visualization_thread.joinable())
        Visualization_thread.join();

    flagExit = false;
    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    if (progressBarThread.joinable())
        progressBarThread.join();

    // progressBarThread.;
}

void Daizy::ShowMB()
{
    std::string tmp = errorMsg;
    errorMsg.clear();
    QMessageBox::critical(this, "Daisi error", tmp.c_str());
};
void stringToQstring(QString& str1, const std::vector<std::string>& str2)
{
    if (!str2.size())
        return;

    for (int i = 0; i < str2.size() - 1; i++)
        str1 = str1 + QString(str2[i].c_str()) + "\n";

    str1 = str1 + QString(str2.back().c_str());
    str1 = str1 + "...\n";
}
void Daizy::setProgress()
{
    QTime time;
    while (flagExit)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        if (errorMsg.size())
        {
            emit showMBSignal();
        }
        if (progress > 0 && progress < 0.01)
            progress = 0.015;
        int progr    = int(progress * 100);
        int progrLoc = int(progressLoc * 100);

        if (progr == 100 || progr == 0)
            time.start();

        int     ms      = time.elapsed();
        int     minuts  = floor(ms / (60 * 1000));
        int     os      = ms - minuts * 60 * 1000;
        int     sec     = floor(os / 1000);
        QString timeStr = QString("time elapsed: ") + QString::number(minuts) + QString(" minuts, ") +
                          QString::number(sec) + QString(" seconds ");
        if (progr != 100)
            QMetaObject::invokeMethod(labelTime, "setText", Qt::QueuedConnection, Q_ARG(QString, timeStr));

        QMetaObject::invokeMethod(progressBar, "setValue", Qt::QueuedConnection, Q_ARG(int, progr));
        QMetaObject::invokeMethod(progressBarLocal, "setValue", Qt::QueuedConnection, Q_ARG(int, progrLoc));
        QString statusL;
        stringToQstring(statusL, status);
        if (statusL != statusLast)
        {
            QMetaObject::invokeMethod(monitor, "setText", Qt::QueuedConnection, Q_ARG(QString, statusL));
            int maxS  = scrollBar->maximum();
            int maxS1 = QMetaObject::invokeMethod(scrollBar, "maximum", Qt::QueuedConnection);

            QMetaObject::invokeMethod(scrollBar, "setValue", Qt::QueuedConnection, Q_ARG(int, 2 * maxS));

            statusLast = statusL;
        }

        /*	if (flagAbort)
                {
                        std::this_thread::sleep_for(5000);
                        for (int i = 0; i < VTKArray.size(); i++)
                        {
                                VTKArray[i]->refresh(0, 1);
                                //QMetaObject::invokeMethod(VTKArray[i], "refresh", Qt::QueuedConnection);
                                //VTKArray[i]->refresh(0);
                        }
                }*/

        //	progressBar->setValue(progress * 100);
    };
};
void Daizy::AbortSimulations()
{
    flagAbort = false;
};
void Daizy::updateGr()
{
    /*	if (flagReset == 1)
            {
                    VTKArray[0]->reset();
                    flagReset = 0;
            }*/
}
void Daizy::stopTimer()
{
    if (timerViz->isActive())
    {
        timerViz->stop();
        //	graphicsWidget = new QVTKWidget();
        //	mainVTK = vtkComponent();
        //	rigthWidgetGrid->addWidget(graphicsWidget, 1, 0);
        clearLayout(rigthWidgetGrid);

        for (int i = 0; i < VTKArray.size(); i++)
            delete VTKArray[i];

        VTKArray.clear();

        graphicsArray.clear();

        if (VTKArray.size() == 0)
        {
            graphicsArray.push_back(new QVTKWidget());
            graphicsArray[0]->setMinimumSize(QSize(0, 0));
            VTKArray.push_back(new vtkComponent());
            rigthWidgetGrid->addWidget(graphicsArray.back(), 0, 0);
        }

        VTKArray[0]->setWidget(graphicsArray[0]);
    }
}