#ifndef PROJECT_H
#define PROJECT_H
//
#include "ModelInterface.h"

#include "AccelModelInterface.h"

#include <stdexcept>
// typedef arma::cx_vec containerFloat;
// typedef arma::vec  containerDouble;
namespace Dproject
{
class project
{
    void LoadFromFile(std::string, std::string version, std::string& errorMsg);

  public:
    int                  problemType;
    int                  precisionType;
    int                  solverType;
    std::string          version;
    ModelInterface*      currentModel;
    AccelModelInterface* accelModel;

    std::string  currentModelFileName;
    std::string  projectName;
    std::wstring boundariesFolder;
    std::wstring resultsFolder;
    std::string  meshParamFileName;
    std::string  projectFolder;
    std::string  currentModelFilePath;

    project(std::string, std::string version, std::string& errorMsg);
    void SaveCurrentModel(std::string, std::string version);
    void SaveProject(std::string version);
    void ModelCreate();
    void ModelLoad(std::string, std::string version, std::string& errorMsg);
    void SaveData(std::string);
    void LoadData(std::string);
    void LoadProject(std::string, std::string version, std::string& errorMsg);
    void ModelDelete();
    project();
    project(int problemTypeIn, int precisionTypeIn, int solverTypeIn, std::string projectNameIn,
            std::string mainFolderIn, std::string version);
};
}
#endif
