#ifndef FilesBrouses_H
#define FilesBrouses_H
#include <QtWidgets/QWidget>
#include <functional>

class MyQTextEdit;
class QLabel;
class QGridLayout;
class QGroupBox;
class QVBoxLayout;
class QPushButton;
class QSignalMapper;
class FilesBrouses : public QWidget
{
    Q_OBJECT
    std::vector<MyQTextEdit*> textEditors;
    std::vector<QLabel*>      textEditorsLabels;
    QGroupBox*                groupBox;
    QVBoxLayout*              vbox;
    QSignalMapper*            FileBrouseMapper;
    std::vector<QPushButton*> pushBrouses;
    std::string               path;
    std::vector<std::string>  brouseExt;
    std::vector<std::function<int(std::string, std::string&)>> setter;

  public:
    QGroupBox* GetPointer();
    FilesBrouses();
    ~FilesBrouses();
    FilesBrouses(const QString& name);
    void Create(const std::vector<std::string>& brouseNames, const std::vector<std::string>& brouseExtIn,
                const std::vector<std::string>& filesNames, const std::string& pathIn,
                std::vector<std::function<int(std::string, std::string&)>>     setterIn);
  public slots:
    void clickedBrouseSlot(const int& i);
    bool GetParameters(std::vector<std::string>& param);
};
#endif
