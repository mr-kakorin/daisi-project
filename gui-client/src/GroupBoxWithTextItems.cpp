#include "GroupBoxWithTextItems.h"
#include "MyTreeItem.h"
#include <QFileDialog>
#include <QRadioButton>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

QGroupBox* GroupBoxWithTextItems::GetPointer()
{
    return groupBox;
};
GroupBoxWithTextItems::GroupBoxWithTextItems(){

};
GroupBoxWithTextItems::GroupBoxWithTextItems(const QString& name)
{
    groupBox = new QGroupBox(name);
    vbox     = new QVBoxLayout();
    groupBox->setLayout(vbox);
};
GroupBoxWithTextItems::~GroupBoxWithTextItems()
{
    for (int i = 0; i < textEditors.size(); i++)
    {
        delete textEditorsLabels[i];
        delete textEditors[i];
    }

    delete vbox;
    delete groupBox;
}

bool GroupBoxWithTextItems::GetParameters(std::vector<std::string>& param)
{
    param.clear();
    for (int j = 0; j < textEditors.size(); j++)
        param.push_back(textEditors[j]->toPlainText().toStdString());

    return true;
};
void GroupBoxWithTextItems::Create(const std::vector<std::string>& textEditorsNames,
                                   const std::vector<std::string>& textEditorsValues)
{
    textEditors.clear();
    textEditorsLabels.clear();
    for (int j = 0; j < textEditorsNames.size(); j++)
    {
        textEditorsLabels.push_back(new QLabel());
        textEditorsLabels.back()->setText(textEditorsNames[j].c_str());
        vbox->addWidget(textEditorsLabels.back());
        textEditors.push_back(new MyQTextEdit());
        textEditors.back()->setText(textEditorsValues[j].c_str());
        textEditors.back()->setMaximumSize(QSize(16777215, 50));
        vbox->addWidget(textEditors.back());
    };
    groupBox->setMaximumHeight(80 * textEditorsNames.size() + 15);
}
