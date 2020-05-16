#include "GroupBoxWithItems.h"
#include "MyTreeItem.h"
#include <QFileDialog>
#include <QMainWindow>
#include <QPalette>
#include <QPixmap>
#include <QRadioButton>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

void GroupBoxWithItems::SetTextBoxEnable(int i, bool flag)
{
    textEditors[i]->setEnabled(flag);
};

void GroupBoxWithItems::SetParameters(const std::vector<double>& parameters)
{
    for (int j = 0; j < textEditors.size(); j++)
        textEditors[j]->setText(QString::number(parameters[j]));
};

void GroupBoxWithItems::radioButtonClick()
{
    listExec();
};

void GroupBoxWithItems::SetRadioButtonsSignals(const std::function<void()>& listExecIn)
{
    listExec = listExecIn;
    for (int i = 0; i < radioButtons.size(); i++)
        connect(radioButtons[i], SIGNAL(clicked()), this, SLOT(radioButtonClick()));
};

void GroupBoxWithItems::SetEditorsDisabled()
{
    for (int j = 0; j < textEditors.size(); j++)
        textEditors[j]->setEnabled(false);
};
QGroupBox* GroupBoxWithItems::GetPointer()
{
    return groupBox;
};
GroupBoxWithItems::GroupBoxWithItems(){

};
GroupBoxWithItems::GroupBoxWithItems(const QString& name)
{
    groupBox = new QGroupBox(name);
    grid     = new QGridLayout();
    groupBox->setLayout(grid);
};
void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<std::string>& textEditorsValues, const std::string& flag)
{
    textEditorsLabels.clear();
    textLabels.clear();
    for (int j = 0; j < textEditorsNames.size(); j++)
    {
        textEditorsLabels.push_back(new QLabel());
        textEditorsLabels.back()->setText(textEditorsNames[j].c_str());
        grid->addWidget(textEditorsLabels.back(), j, 0);
        textLabels.push_back(new QLabel());
        textLabels.back()->setText(textEditorsValues[j].c_str());
        grid->addWidget(textLabels.back(), j, 1);
    };
    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 30);
};

void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<double>& textEditorsValues, const std::string& flag)
{
    textEditorsLabels.clear();
    textLabels.clear();
    for (int j = 0; j < textEditorsNames.size(); j++)
    {
        textEditorsLabels.push_back(new QLabel());
        textEditorsLabels.back()->setText(textEditorsNames[j].c_str());
        grid->addWidget(textEditorsLabels.back(), j, 0);
        textLabels.push_back(new QLabel());
        textLabels.back()->setText(QString::number(textEditorsValues[j]));
        grid->addWidget(textLabels.back(), j, 1);
    };
    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 30);
};

void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<double>&      textEditorsValues)
{
    textEditorsLabels.clear();
    textEditors.clear();
    for (int j = 0; j < textEditorsNames.size(); j++)
    {
        textEditorsLabels.push_back(new QLabel());
        textEditorsLabels.back()->setText(textEditorsNames[j].c_str());
        grid->addWidget(textEditorsLabels.back(), j, 0);
        textEditors.push_back(new MyQTextEdit());
        textEditors.back()->setMaximumSize(QSize(16777215, 28));
        textEditors.back()->setText(QString::number(textEditorsValues[j]));
        connect(textEditors.back(), SIGNAL(textChanged()), this, SLOT(checkInput()));
        grid->addWidget(textEditors.back(), j, 1);
    };
    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 30);
}

void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<int>&         textEditorsValues)
{
    textEditorsLabels.clear();
    textEditors.clear();
    for (int j = 0; j < textEditorsNames.size(); j++)
    {
        textEditorsLabels.push_back(new QLabel());
        textEditorsLabels.back()->setText(textEditorsNames[j].c_str());
        grid->addWidget(textEditorsLabels.back(), j, 0);
        textEditors.push_back(new MyQTextEdit());
        textEditors.back()->setMaximumSize(QSize(16777215, 28));
        textEditors.back()->setText(QString::number(textEditorsValues[j]));
        connect(textEditors.back(), SIGNAL(textChanged()), this, SLOT(checkInput()));
        grid->addWidget(textEditors.back(), j, 1);
    };
    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 30);
}

void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<double>&      textEditorsValues,
                               const std::vector<std::string>& radioButtonsNames, double radioButtonsValue,
                               std::vector<int> flags)
{
    if (!flags.size())
    {
        flags.resize(radioButtonsNames.size());

        for (int i   = 0; i < radioButtonsNames.size(); i++)
            flags[i] = 1;
    }

    Create(textEditorsNames, textEditorsValues);
    int icurr = textEditorsNames.size();
    radioButtons.clear();
    for (int i = 0; i < radioButtonsNames.size(); i++)
    {
        radioButtons.push_back(new QRadioButton(radioButtonsNames[i].c_str()));
        grid->addWidget(radioButtons.back(), icurr + i, 0);
        radioButtons.back()->setChecked(false);
        radioButtons.back()->setEnabled(flags[i]);
    }
    if (radioButtons.size())
        radioButtons[int(radioButtonsValue)]->setChecked(true);

    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 35 * radioButtonsNames.size() + 30);
}
void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<double>&      textEditorsValues,
                               const std::vector<std::string>& radioButtonsNames, double radioButtonsValue,
                               const std::vector<std::string>& checkBoxesNames,
                               const std::vector<double>& checkBoxesValues, std::vector<int> flags)
{

    Create(textEditorsNames, textEditorsValues, radioButtonsNames, radioButtonsValue, flags);

    int icurr = textEditorsNames.size() + radioButtonsNames.size();
    checkBoxes.clear();

    for (int i = 0; i < checkBoxesNames.size(); i++)
    {
        checkBoxes.push_back(new QCheckBox(checkBoxesNames[i].c_str()));
        checkBoxes.back()->setChecked(int(checkBoxesValues[i]));
        grid->addWidget(checkBoxes.back(), icurr + i, 0);
    }
    groupBox->setMaximumHeight(45 * textEditorsNames.size() + 25 * checkBoxesNames.size() +
                               35 * radioButtonsNames.size() + 30);
}
void GroupBoxWithItems::Create(const std::vector<std::string>& textEditorsNames,
                               const std::vector<std::string>& radioButtonsNames,
                               const std::vector<std::string>& checkBoxesNames, const std::vector<double>& parameters,
                               std::vector<int> flags)
{
    Create(textEditorsNames, std::vector<double>(parameters.begin(), parameters.begin() + textEditorsNames.size()),
           radioButtonsNames, parameters[textEditorsNames.size()], checkBoxesNames,
           std::vector<double>(parameters.begin() + textEditorsNames.size() + 1, parameters.end()), flags);
}

GroupBoxWithItems::~GroupBoxWithItems()
{
    for (int i = 0; i < textEditors.size(); i++)
    {
        delete textEditors[i];
        delete textEditorsLabels[i];
    }

    for (int i = 0; i < checkBoxes.size(); i++)
        delete checkBoxes[i];

    for (int i = 0; i < radioButtons.size(); i++)
        delete radioButtons[i];

    for (int i = 0; i < textLabels.size(); i++)
        delete textLabels[i];

    delete grid;
    delete groupBox;
};
void GroupBoxWithItems::checkInput()
{
    bool ok = true;
    int  s  = textEditors.size();

    for (int j = 0; j < s; j++)
        textEditors[j]->toDoubleMy(&ok);

    emit textChanged();
};
bool GroupBoxWithItems::GetParameters(std::vector<double>& tmp3)
{
    bool ok = true;
    int  s  = textEditors.size();
    tmp3.clear();

    for (int j = 0; j < s; j++)
    {
        bool locok;
        tmp3.push_back(textEditors[j]->toDoubleMy(&locok));
        if (!locok)
        {
            ok = false;
        };
    }

    for (int i = 0; i < radioButtons.size(); i++)
    {
        if (radioButtons[i]->isChecked())
        {
            tmp3.push_back(i);
            break;
        }
    };

    for (int i = 0; i < checkBoxes.size(); i++)
        tmp3.push_back(checkBoxes[i]->isChecked());

    return ok;
    /*if (!ok)
    {
            QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }*/
};