#ifndef AccelModelInterface_H
#define AccelModelInterface_H

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class SimulationDataAccel;
class AccelModelInterface
{
  public:
    virtual void SetDataFolder(const std::string& dataFolderIn) = 0;
    virtual void SetAccelSectionParameters(int                        currentAccelSection,
                                           const std::vector<double>& params) = 0;
    virtual std::vector<std::string> GetMainAccelParameterNamesCalculated()   = 0;
    virtual int SetSomeParametersFromFile(int n, const std::string& filename,
                                          std::string& errorMessage) = 0;
    virtual int SetSomeSolverFileName(const std::string& solverName, int fileNumber, int io,
                                      const std::string& filename, std::string& errorMessage) = 0;
    virtual int SetSolverAllParameters(const std::string&              solverName,
                                       const std::vector<std::string>& Inputfilenames,
                                       const std::vector<std::string>& Outputfilenames,
                                       const std::vector<double>&      in,
                                       const std::vector<double>&      inF,
                                       std::string&                    errorMessage) = 0;
    virtual void GetSomeSolverFileName(const std::string&        solverName,
                                       std::vector<std::string>& input,
                                       std::vector<std::string>& output,
                                       std::vector<std::string>& filesInDescr,
                                       std::vector<std::string>& filesOutDescr) = 0;
    virtual void Solve(const std::string& solverName, double& progress, bool& flagAbort,
                       std::string& errorMsg, std::vector<std::string>& status) = 0;
    virtual std::vector<std::string> GetSomeFileName() = 0;
    virtual void GetBrouseFlags(std::vector<std::string>& brouse,
                                std::vector<std::string>& names) = 0;
    virtual size_t GetNumberOfAccelFlows()                       = 0;
    virtual void   AddFlow()                                     = 0;
    virtual void GetParametersAccelFlow(std::vector<std::string>& keys, std::vector<double>& p,
                                        int flowNumber) = 0;
    virtual void SetParametersAccelFlow(const std::vector<double>& params, int flowNumber) = 0;
    virtual void GetMainAccelParameters(std::vector<std::string>& keysO,
                                        std::vector<double>&      p) = 0;
    virtual void GetMainAccelParametersFlags(std::vector<std::string>& keysO,
                                             std::vector<double>&      p) = 0;

    virtual std::vector<double> GetMainAccelParameterCalculated()      = 0;
    virtual void SetMainAccelParameters(const std::vector<double>& in) = 0;

    virtual void SetSolverParameters(const std::string&         solverName,
                                     const std::vector<double>& in) = 0;
    virtual void GetSolverParameters(const std::string& solverName, std::vector<std::string>& keysO,
                                     std::vector<double>& p) = 0;
    virtual void GetSolverParametersFlags(const std::string&                     solverName,
                                          std::vector<std::vector<std::string>>& keysO,
                                          std::vector<double>&                   p) = 0;
    virtual void SetSolverParametersFlags(const std::string&         solverName,
                                          const std::vector<double>& in) = 0;

    virtual void GetSectionParameters(std::vector<std::vector<std::string>>& keysO,
                                      std::vector<std::vector<double>>&      p) = 0;

    virtual void GetSolverFileProps(const std::string& solverName, int fileNumber, std::string& ext,
                                    std::string& prop)                             = 0;
    virtual std::vector<std::shared_ptr<SimulationDataAccel>>& GetSimulationData() = 0;
    //	virtual std::vector<SimulationDataAccel*>& GetSimulationData() = 0;
    virtual void GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                            std::vector<std::string>&         names) = 0;
    virtual std::vector<std::string> GetSolversNames()                               = 0;
    virtual std::vector<std::string> GetnamesOfSequences()                           = 0;
    virtual std::vector<std::string> GetcalcParametersNames()                        = 0;
    virtual void GetSequencesOfParameters(int seq, std::vector<double>& sequence)     = 0;
    virtual void GetSequencesOfParametersCalc(int seq, std::vector<double>& sequence) = 0;
    virtual void SaveData(std::string fileName) = 0;
    virtual void LoadData(std::string fileName) = 0;
};
#endif