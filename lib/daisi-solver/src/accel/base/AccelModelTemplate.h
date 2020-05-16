#ifndef AccelModelTemplate_H
#define AccelModelTemplate_H

#include "AccelModelInterface.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

//#include "../base/AccelFlow.h"
template <class DeviceType, class SolverType>
class AccelModelTemplate : public AccelModelInterface
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  protected:
    std::shared_ptr<DeviceType>                       device;
    std::shared_ptr<SolverType>                       solver;
    std::string                                       dataFolder;
    std::vector<std::shared_ptr<SimulationDataAccel>> outputData;

  public:
    void SetDataFolder(const std::string& dataFolderIn)
    {
        dataFolder = dataFolderIn;
    };
    void SaveData(std::string fileName);
    void LoadData(std::string fileName);
    AccelModelTemplate(const std::string& dataFolderIn);
    AccelModelTemplate(){};
    //~AccelModelTemplate();
    std::vector<std::shared_ptr<SimulationDataAccel>>& GetSimulationData();
    /*std::vector<SimulationDataAccel*>& GetSimulationData()
    {
            return dd;
    };*/
    void SetAccelSectionParameters(int currentAccelSection, const std::vector<double>& params);
    void GetSectionParameters(std::vector<std::vector<std::string>>& keysO,
                              std::vector<std::vector<double>>&      p);
    size_t GetNumberOfAccelFlows();
    int SetSomeParametersFromFile(int n, const std::string& filename, std::string& errorMessage);
    int SetSomeSolverFileName(const std::string& solverName, int fileNumber, int io,
                              const std::string& filename, std::string& errorMessage);
    void GetSomeSolverFileName(const std::string& solverName, std::vector<std::string>& input,
                               std::vector<std::string>& output,
                               std::vector<std::string>& filesInDescr,
                               std::vector<std::string>& filesOutDescr);
    void GetSolverFileProps(const std::string& solverName, int fileNumber, std::string& ext,
                            std::string& prop);
    std::vector<std::string> GetSomeFileName();
    void                     AddFlow();
    void GetParametersAccelFlow(std::vector<std::string>& keys, std::vector<double>& p,
                                int flowNumber);
    void SetParametersAccelFlow(const std::vector<double>& params, int flowNumber);
    void GetMainAccelParameters(std::vector<std::string>& keysO, std::vector<double>& p);
    void GetMainAccelParametersFlags(std::vector<std::string>& keysO, std::vector<double>& p);

    std::vector<double> GetMainAccelParameterCalculated();
    void SetMainAccelParameters(const std::vector<double>& in);
    void SetSolverParameters(const std::string& solverName, const std::vector<double>& in);
    void GetSolverParameters(const std::string& solverName, std::vector<std::string>& keysO,
                             std::vector<double>& p);

    void SetSolverParametersFlags(const std::string& solverName, const std::vector<double>& in);
    void GetSolverParametersFlags(const std::string&                     solverName,
                                  std::vector<std::vector<std::string>>& keysO,
                                  std::vector<double>&                   p);

    std::vector<std::string> GetSolversNames();
    void GetBrouseFlags(std::vector<std::string>& brouse, std::vector<std::string>& names);
    std::vector<std::string> GetnamesOfSequences();
    void GetSequencesOfParameters(int seq, std::vector<double>& sequence);
    void GetSequencesOfParametersCalc(int seq, std::vector<double>& sequence);
    std::vector<std::string> GetcalcParametersNames();
    std::vector<std::string> GetMainAccelParameterNamesCalculated();
    int SetSolverAllParameters(const std::string&              solverName,
                               const std::vector<std::string>& Inputfilenames,
                               const std::vector<std::string>& Outputfilenames,
                               const std::vector<double>& in, const std::vector<double>& inF,
                               std::string& errorMessage);
};

#endif