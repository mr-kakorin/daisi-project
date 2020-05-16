#include "FilesBrouses.h"
#include "MyTreeItem.h"
#include <QSignalMapper>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

QGroupBox* FilesBrouses::GetPointer()
{
    return groupBox;
};
FilesBrouses::FilesBrouses(){

};
FilesBrouses::~FilesBrouses()
{
    for (int i = 0; i < textEditors.size(); i++)
    {
        delete textEditors[i];
        delete textEditorsLabels[i];
        delete pushBrouses[i];
    }

    delete vbox;
    delete groupBox;
};
bool FilesBrouses::GetParameters(std::vector<std::string>& param)
{
    param.clear();
    for (int j = 0; j < textEditors.size(); j++)
        param.push_back(textEditors[j]->toPlainText().toStdString());

    return true;
};
FilesBrouses::FilesBrouses(const QString& name)
{
    groupBox         = new QGroupBox(name);
    vbox             = new QVBoxLayout;
    FileBrouseMapper = new QSignalMapper();
    groupBox->setLayout(vbox);
};
void FilesBrouses::Create(const std::vector<std::string>& brouseNames, const std::vector<std::string>& brouseExtIn,
                          const std::vector<std::string>& filesNames, const std::string& pathIn,
                          std::vector<std::function<int(std::string, std::string&)>>     setterIn)
{
    setter    = setterIn;
    path      = pathIn;
    brouseExt = brouseExtIn;
    groupBox->setMaximumHeight(brouseNames.size() * 110);

    textEditors.resize(brouseNames.size());
    textEditorsLabels.resize(brouseNames.size());
    pushBrouses.resize(brouseNames.size());
    for (int i = 0; i < brouseNames.size(); i++)
    {
        textEditorsLabels[i] = new QLabel(brouseNames[i].c_str());

        vbox->addWidget(textEditorsLabels[i]);

        textEditors[i] = new MyQTextEdit();
        textEditors[i]->setMaximumSize(QSize(16777215, 50));
        textEditors[i]->setText(filesNames[i].c_str());

        vbox->addWidget(textEditors[i]);

        pushBrouses[i] = new QPushButton("Browse");
        pushBrouses[i]->setMaximumSize(QSize(60, 40));
        connect(pushBrouses[i], SIGNAL(clicked()), FileBrouseMapper, SLOT(map()));
        FileBrouseMapper->setMapping(pushBrouses[i], i);
        vbox->addWidget(pushBrouses[i]);
    }

    connect(FileBrouseMapper, SIGNAL(mapped(const int&)), this, SLOT(clickedBrouseSlot(const int&)));
};

void FilesBrouses::clickedBrouseSlot(const int& i)
{
    std::string errorMessage;

    QString fileName = QFileDialog::getOpenFileName(this, "Open File", path.c_str(), brouseExt[i].c_str());

    if (!fileName.size())
        return;

    textEditors[i]->setText(fileName);

    int succes;
    //	if (fileName.size())
    //		succes = (*currentProject)->accelModel->SetSomeParametersFromFile(i, fileName.toStdString(),
    // errorMessage);
    succes = setter[i](fileName.toStdString(), errorMessage);

    if (!succes)
    {
        QMessageBox messageBox;
        messageBox.critical(0, "Error", errorMessage.c_str());
        messageBox.setFixedSize(500, 200);
    }
    //	else
    //	showAccelParameters(0);
};