#ifndef GroupBoxWithTextItems_H
#define GroupBoxWithTextItems_H
#include <QtWidgets/QWidget>
class QGroupBox;
class QVBoxLayout;
class MyQTextEdit;
class QLabel;
class GroupBoxWithTextItems : public QWidget
{
    Q_OBJECT
    std::vector<MyQTextEdit*> textEditors;
    std::vector<QLabel*>      textEditorsLabels;
    QGroupBox*                groupBox;
    QVBoxLayout*              vbox;

  public:
    QGroupBox* GetPointer();
    GroupBoxWithTextItems();
    GroupBoxWithTextItems(const QString& name);
    ~GroupBoxWithTextItems();
    bool GetParameters(std::vector<std::string>& param);
    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<std::string>& textEditorsValues);
};
#endif
