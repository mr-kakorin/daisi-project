#ifndef GroupBoxWithItems_H
#define GroupBoxWithItems_H
#include <QtWidgets/QWidget>
#include <functional>
class MyQTextEdit;
class QLabel;
class QCheckBox;
class QRadioButton;
class QGridLayout;
class QGroupBox;
class GroupBoxWithItems : public QWidget
{
    Q_OBJECT
    std::vector<MyQTextEdit*> textEditors;
    std::vector<QLabel*>      textEditorsLabels;
    std::vector<QLabel*>      textLabels;

    std::vector<QCheckBox*>    checkBoxes;
    std::vector<QRadioButton*> radioButtons;
    QGroupBox*                 groupBox;
    QGridLayout*               grid;
    std::function<void()>      listExec;

  public:
    void SetTextBoxEnable(int i, bool flag);
    QGroupBox* GetPointer();
    GroupBoxWithItems();
    GroupBoxWithItems(const QString& name);
    void SetParameters(const std::vector<double>& parameters);
    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<double>& textEditorsValues);
    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<int>& textEditorsValues);

    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<double>& textEditorsValues,
                const std::string& flag);

    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<std::string>& textEditorsValues,
                const std::string& flag);

    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<double>& textEditorsValues,
                const std::vector<std::string>& radioButtonsNames, double radioButtonsValue,
                std::vector<int> flags = {});

    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<double>& textEditorsValues,
                const std::vector<std::string>& radioButtonsNames, double radioButtonsValue,
                const std::vector<std::string>& checkBoxesNames, const std::vector<double>& checkBoxesValues,
                std::vector<int> flags = {});

    void Create(const std::vector<std::string>& textEditorsNames, const std::vector<std::string>& radioButtonsNames,
                const std::vector<std::string>& checkBoxesNames, const std::vector<double>& parameters,
                std::vector<int> flags = {});

    ~GroupBoxWithItems();
    void SetEditorsDisabled();
    bool GetParameters(std::vector<double>& param);
    void SetRadioButtonsSignals(const std::function<void()>& listExecIn);
  public slots:
    void checkInput();
    void radioButtonClick();
  signals:
    void textChanged();
};
#endif
