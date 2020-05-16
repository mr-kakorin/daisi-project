#ifndef AccelSolverBase_H
#define AccelSolverBase_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <map>
class myunsorted_map;
class AccelSolverBase
{

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const;
    template <class Archive> void load(Archive& ar, const unsigned int);

  protected:
    std::vector<std::string> solversNames;
    std::map<std::string, std::vector<std::string>> inputFileNames;
    std::map<std::string, std::vector<std::string>> inputFileExt;
    std::map<std::string, std::vector<std::string>> inputFileProps;

    std::map<std::string, std::vector<std::string>> outpuFileNames;
    std::map<std::string, std::vector<std::string>> outpuFilePaths;

    std::map<std::string, std::vector<std::string>> outpuFileDescr;
    std::map<std::string, std::vector<std::string>> inputFileDescr;

    std::map<std::string, std::shared_ptr<myunsorted_map>>       solverParameters;
    std::map<std::string, std::vector<std::vector<std::string>>> solverParametersFlagsNames;
    std::map<std::string, std::vector<double>>                   solverParametersFlags;

    void BaseInit(const std::string& dataFolderIn);

  public:
    void ChangeSolver(const std::string&                           solverName,
                      const std::vector<std::vector<std::string>>& description,
                      std::vector<std::string>&                    solverParametersNames,
                      std::vector<std::vector<std::string>>&       solverParametersNamesFlags);
    void AddSolver(const std::string&                           solverName,
                   const std::vector<std::vector<std::string>>& description,
                   std::vector<std::string>&                    solverParametersNames,
                   std::vector<std::vector<std::string>>&       solverParametersNamesFlags);
    void SetFolder(const std::string& dataFolderIn);
    void checkInputParameters(std::string& dataFolder, std::string& errorMessage,
                              const std::string& solverName);
    AccelSolverBase(const std::string& dataFolderIn, int type);
    AccelSolverBase(){};
    int SetSomeSolverFileName(const std::string& solverName, int fileNumber, int io,
                              const std::string& filename, std::string& errorMessage);
    void GetSomeSolverFileName(const std::string& solverName, std::vector<std::string>& input,
                               std::vector<std::string>& output,
                               std::vector<std::string>& filesInDescr,
                               std::vector<std::string>& filesOutDescr);
    void GetSolverFileProps(const std::string& solverName, int fileNumber, std::string& ext,
                            std::string& prop);
    void SetSolverParameters(const std::string& solverName, const std::vector<double>& in);
    void GetSolverParameters(const std::string& solverName, std::vector<std::string>& keysO,
                             std::vector<double>& p);
    void SetSolverParametersFlags(const std::string& solverName, const std::vector<double>& in);
    void GetSolverParametersFlags(const std::string&                     solverName,
                                  std::vector<std::vector<std::string>>& keysO,
                                  std::vector<double>&                   p);
    std::vector<std::string> GetSolversNames();
};
#endif
