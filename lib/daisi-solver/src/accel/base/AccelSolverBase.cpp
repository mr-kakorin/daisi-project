#include "AccelSolverBase.h"
#include "Results.h"
#include "Tools.h"

template void
AccelSolverBase::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                            const unsigned int file_version);
template void
AccelSolverBase::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                            const unsigned int file_version);

template <class Archive> void AccelSolverBase::save(Archive& ar, const unsigned int) const
{
    ar& inputFileNames;
    ar& outpuFileNames;
    ar& inputFileExt;
    ar& inputFileProps;
    ar& solverParameters;
    ar& solversNames;
    ar& inputFileDescr;
    ar& outpuFileDescr;
    ar& solverParametersFlags;
    ar& solverParametersFlagsNames;
};
template <class Archive> void AccelSolverBase::load(Archive& ar, const unsigned int)
{
    ar& inputFileNames;
    ar& outpuFileNames;
    ar& inputFileExt;
    ar& inputFileProps;
    ar& solverParameters;
    ar& solversNames;
    ar& inputFileDescr;
    ar& outpuFileDescr;
    ar& solverParametersFlags;
    ar& solverParametersFlagsNames;
};
void AccelSolverBase::ChangeSolver(
    const std::string& solverName, const std::vector<std::vector<std::string>>& description,
    std::vector<std::string>&              solverParametersNames,
    std::vector<std::vector<std::string>>& solverParametersNamesFlags)
{
    inputFileNames[solverName] = description[0];
    outpuFileNames[solverName] = description[1];
    inputFileExt[solverName]   = description[2];
    inputFileProps[solverName] = description[3];
    outpuFileDescr[solverName] = description[4];
    inputFileDescr[solverName] = description[5];

    solverParameters[solverName] = std::shared_ptr<myunsorted_map>(new myunsorted_map());

    solverParameters[solverName]->clear();

    for (int i = 0; i < solverParametersNames.size(); i++)
        solverParameters[solverName]->insert(solverParametersNames[i], 0.0);

    solverParametersFlagsNames[solverName] = solverParametersNamesFlags;

    solverParametersFlags[solverName] = std::vector<double>{};

    solverParametersFlags[solverName].clear();
    for (int i = 0; i < solverParametersNamesFlags.size(); i++)
        solverParametersFlags[solverName].push_back(1.0);
}
void AccelSolverBase::AddSolver(const std::string&                           solverName,
                                const std::vector<std::vector<std::string>>& description,
                                std::vector<std::string>&                    solverParametersNames,
                                std::vector<std::vector<std::string>>& solverParametersNamesFlags)
{
    solversNames.push_back(solverName);
    inputFileNames.insert(std::make_pair(solverName, description[0]));
    outpuFileNames.insert(std::make_pair(solverName, description[1]));
    inputFileExt.insert(std::make_pair(solverName, description[2]));
    inputFileProps.insert(std::make_pair(solverName, description[3]));
    outpuFileDescr.insert(std::make_pair(solverName, description[4]));
    inputFileDescr.insert(std::make_pair(solverName, description[5]));

    solverParameters.insert(
        std::make_pair(solverName, std::shared_ptr<myunsorted_map>(new myunsorted_map())));

    for (int i = 0; i < solverParametersNames.size(); i++)
        solverParameters[solverName]->insert(solverParametersNames[i], 0.0);

    solverParametersFlagsNames.insert(std::make_pair(solverName, solverParametersNamesFlags));

    solverParametersFlags.insert(std::make_pair(solverName, std::vector<double>{}));

    for (int i = 0; i < solverParametersNamesFlags.size(); i++)
        solverParametersFlags[solverName].push_back(1.0);
};

void AccelSolverBase::checkInputParameters(std::string& dataFolder, std::string& errorMessage,
                                           const std::string& solverName)
{

    SetFolder(dataFolder);

    FILE* fp;

    for (int j = 0; j < inputFileNames[solverName].size(); j++)
    {
        // if(inputFileNames[i][j])
        fp = fopen(inputFileNames[solverName][j].c_str(), "r");
        if (!fp)
        {
            errorMessage =
                errorMessage + "Unable to open input file. " + inputFileNames[solverName][j];
            return;
        }
        fclose(fp);
        fp = nullptr;
    }

    for (int j = 0; j < outpuFilePaths[solverName].size(); j++)
    {
        fp = fopen(outpuFilePaths[solverName][j].c_str(), "w");
        if (!fp)
        {
            errorMessage =
                errorMessage + "Unable to open output file. " + outpuFilePaths[solverName][j];
            return;
        };
        fclose(fp);
        fp = nullptr;
    }
};

AccelSolverBase::AccelSolverBase(const std::string& dataFolderIn, int type){
    /*inputFileNames = flagStringsSolver::SynchrotronInputFiles;
    outpuFilePaths = flagStringsSolver::SynchrotronOutputFiles;
    inputFileExt = flagStringsSolver::SynchrotronInputFilesExt;
    inputFileProps = flagStringsSolver::SynchrotronInputFilesProps;

    for (int i = 0; i<4; i++)
            solverParameters.push_back(std::shared_ptr<myunsorted_map>(new myunsorted_map()));

    solversNames = flagStringsSolver::SynchrotronSolversNames;
    for (int j = 0; j < solverParameters.size(); j++)
    {
            for (int i = 0; i < flagStringsSolver::SynchrotronSolverParameteresNames[j].size(); i++)
                    solverParameters[j]->insert(flagStringsSolver::SynchrotronSolverParameteresNames[j][i],
    0.0);
    }*/
};

void AccelSolverBase::SetFolder(const std::string& dataFolderIn)
{
    outpuFilePaths = outpuFileNames;
    for (auto it = outpuFilePaths.begin(); it != outpuFilePaths.end(); ++it)
    {
        for (int j = 0; j < it->second.size(); j++)
        {
            if (it->second[j].size() != 0)
                it->second[j] = dataFolderIn + it->second[j];
        }
    }
}

int AccelSolverBase::SetSomeSolverFileName(const std::string& solverName, int fileNumber, int io,
                                           const std::string& filename, std::string& errorMessage)
{
    if (io)
        inputFileNames[solverName][fileNumber] = filename;
    else
        outpuFileNames[solverName][fileNumber] = filename;

    return 1;
};
void AccelSolverBase::GetSomeSolverFileName(const std::string&        solverName,
                                            std::vector<std::string>& input,
                                            std::vector<std::string>& output,
                                            std::vector<std::string>& filesInDescr,
                                            std::vector<std::string>& filesOutDescr)
{
    input         = inputFileNames[solverName];
    output        = outpuFileNames[solverName];
    filesInDescr  = inputFileDescr[solverName];
    filesOutDescr = outpuFileDescr[solverName];
};
void AccelSolverBase::SetSolverParameters(const std::string&         solverName,
                                          const std::vector<double>& in)
{
    solverParameters[solverName]->SetValues(in);
};
void AccelSolverBase::SetSolverParametersFlags(const std::string&         solverName,
                                               const std::vector<double>& in)
{
    solverParametersFlags[solverName] = in;
};
void AccelSolverBase::GetSolverParameters(const std::string&        solverName,
                                          std::vector<std::string>& keysO, std::vector<double>& p)
{
    keysO = solverParameters[solverName]->GetKeys();
    p     = solverParameters[solverName]->GetValues();
};
void AccelSolverBase::GetSolverParametersFlags(const std::string&                     solverName,
                                               std::vector<std::vector<std::string>>& keysO,
                                               std::vector<double>&                   p)
{
    keysO = solverParametersFlagsNames[solverName];
    p     = solverParametersFlags[solverName];
};
std::vector<std::string> AccelSolverBase::GetSolversNames()
{
    return solversNames;
};

void AccelSolverBase::GetSolverFileProps(const std::string& solverName, int fileNumber,
                                         std::string& ext, std::string& prop)
{
    ext  = inputFileExt[solverName][fileNumber];
    prop = inputFileProps[solverName][fileNumber];
    // solver->GetSomeSolverFileName( solverName, ext, prop);
};
void AccelSolverBase::BaseInit(const std::string& dataFolderIn)
{
    SetFolder(dataFolderIn);
};