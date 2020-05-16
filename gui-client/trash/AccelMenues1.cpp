#include "DaisiClient.h"
#include "FlagStringsD.h"
std::string GetFileName(const std::string& fullPath)
{
    std::string result = "";
    int         i      = fullPath.size() - 1;
    if (i == -1)
        return result;
    while (1)
    {
        if ((fullPath[i] == '/' && fullPath[i - 1] == '/') || (fullPath[i] == '\\'))
        {
            result = fullPath.substr(i + 1, fullPath.size() - 1);
            return result;
        };
        i--;
        if (i == 0)
            break;
    };
    return result;
};

void Daizy::SetAccelFlowParameters()
{
    bool                ok = true;
    std::vector<double> params;
    for (int i = 0; i < MiddleWidgetlabelsVector.size(); i++)
    {
        params.push_back(MiddleWidgetTextEditVector[i]->toDoubleMy(&ok));
    };
    currentProject->accelModel->SetParametersAccelFlow(params, currentFlow);
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
};

void Daizy::SetAccelSolverParameters()
{
    bool ok = true;

    std::vector<double> params;
    for (int i = 0; i < MiddleWidgetlabelsVector.size(); i++)
    {
        params.push_back(MiddleWidgetTextEditVector[i]->toDoubleMy(&ok));
    };
    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    currentProject->accelModel->SetSolverParameters(currentsolver, params);
};

void Daizy::SetOutputSolverFiles()
{
    for (int j = 0; j < MiddleWidgetTextEditVector1.size(); j++)
        currentProject->accelModel->SetSomeSolverFileName(
            currentProject->projectFolder + "//dataFiles//" +
                MiddleWidgetTextEditVector1[j]->toPlainText().toStdString(),
            currentsolver, j, 0);
}

void Daizy::clickedBrouseSlot(const int& i)
{
    /*MyQTextEdit* editor;
    QString ext;
    std::string path;

    //std::function<int(std::string, std::string)> setter;
    std::string errorMessage;
    path = currentProject->projectFolder;
    std::string flag;
    //std::function<void(AccelModelInterface&, std::string)> setter1 = &AccelModelInterface::SetSequenceFile;
    for (int i = 0; i < flagStrings::brouseFlags.size(); i++)
    {
            if (flagStrings::brouseFlags[i] == text)
            {
                    switch (i)
                    {
                    case 0:
                            editor = MiddleWidgetTextEdit1;
                            ext = flagStrings::FilesExtensions[0];
                            //	setter = std::bind(&AccelModelInterface::SetOpticsFile, currentProject->accelModel,
    std::placeholders::_1, std::placeholders::_2); path = currentProject->projectFolder; break; case 1: editor =
    MiddleWidgetTextEdit2; ext = flagStrings::FilesExtensions[1]; break;
                    };
                    break;
            };
    };*/

    std::string path = currentProject->projectFolder;
    std::string errorMessage;

    QString fileName = QFileDialog::getOpenFileName(this, "Open File", path.c_str(), brouseExt[i].c_str());
    MiddleWidgetTextEditVectorFiles[i]->setText(fileName);

    int succes;
    if (fileName.size())
        succes = currentProject->accelModel->SetSomeParametersFromFile(i, fileName.toStdString(), errorMessage);
    // setter(fileName.toStdString(), errorMessage);

    if (!succes)
    {
        QMessageBox messageBox;
        messageBox.critical(0, "Error", errorMessage.c_str());
        messageBox.setFixedSize(500, 200);
    }
    else
        showAccelParameters(0);
};