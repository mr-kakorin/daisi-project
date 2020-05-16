#pragma once

#ifndef WIN32
#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>
#endif

#ifdef WIN32
#include <daisi-solver/include/Results.h>
#include <daisi-solver/include/project.h>
#endif

#include "MiddleWidget.h"
#include <QFileDialog>
#include <QMainWindow>
#include <QRadioButton>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

class Controller : public QWidget
{
    Q_OBJECT
    MiddleWidget* middleWidget;
    QTreeWidget*  leftWidgetTree;

  public:
    Dproject::project* currentProject;
    void ShowProjectTreeAccel(bool flagExpand = false);
  public slots:
    void LoadProjectEvent();
    void CreateNewProjectEvent();
};