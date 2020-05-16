#include "AccelModelTemplate.h"

#include "../base/AccelFlow.h"
#include "../base/LinacFlow.h"
#include "Results.h"

#ifdef NUCL

#include "../Synchrotron/SynchrotronDevice.h"
#include "../Synchrotron/SynchrotronSolver.h"

template class AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>;

template void AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);

template void
AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template void
AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

#endif

#include "../RFQ/RFQDevice.h"
#include "../RFQ/RFQSolver.h"

template class AccelModelTemplate<RFQDevice, RFQSolver>;

template void AccelModelTemplate<RFQDevice, RFQSolver>::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void AccelModelTemplate<RFQDevice, RFQSolver>::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version);

template void AccelModelTemplate<RFQDevice, RFQSolver>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template void AccelModelTemplate<RFQDevice, RFQSolver>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class DeviceType, class SolverType>
template <class Archive>
void AccelModelTemplate<DeviceType, SolverType>::save(Archive& ar, const unsigned int) const
{
    ar& device;
    ar& solver;
}
template <class DeviceType, class SolverType>
template <class Archive>
void AccelModelTemplate<DeviceType, SolverType>::load(Archive& ar, const unsigned int)
{
    ar& device;
    ar& solver;
}
template <class DeviceType, class SolverType>
int AccelModelTemplate<DeviceType, SolverType>::SetSolverAllParameters(
    const std::string& solverName, const std::vector<std::string>& Inputfilenames,
    const std::vector<std::string>& Outputfilenames, const std::vector<double>& in,
    const std::vector<double>& inF, std::string& errorMessage)
{
    for (int i = 0; i < Inputfilenames.size(); i++)
        SetSomeSolverFileName(solverName, i, 1, Inputfilenames[i], errorMessage);

    for (int i = 0; i < Outputfilenames.size(); i++)
        SetSomeSolverFileName(solverName, i, 0, Outputfilenames[i], errorMessage);

    SetSolverParameters(solverName, in);
    SetSolverParametersFlags(solverName, inF);

    return 1;
};
template <class DeviceType, class SolverType>
std::vector<std::string>
AccelModelTemplate<DeviceType, SolverType>::GetMainAccelParameterNamesCalculated()
{
    return device->GetcalcMainParameters();
};

template <class DeviceType, class SolverType>
std::vector<std::string> AccelModelTemplate<DeviceType, SolverType>::GetcalcParametersNames()
{
    return device->GetcalcParametersNames();
};

template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSequencesOfParameters(
    int seq, std::vector<double>& sequence)
{
    device->GetSequencesOfParameters(seq, sequence);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSequencesOfParametersCalc(
    int seq, std::vector<double>& sequence)
{
    device->GetSequencesOfParametersCalc(seq, sequence);
};
template <class DeviceType, class SolverType>
std::vector<std::string> AccelModelTemplate<DeviceType, SolverType>::GetnamesOfSequences()
{
    return device->GetnamesOfSequences();
};

template <class DeviceType, class SolverType>
std::vector<std::shared_ptr<SimulationDataAccel>>&
AccelModelTemplate<DeviceType, SolverType>::GetSimulationData()
{
    return outputData;
};
template <class DeviceType, class SolverType>
size_t AccelModelTemplate<DeviceType, SolverType>::GetNumberOfAccelFlows()
{
    return device->GetLinacFlows().size();
};
template <class DeviceType, class SolverType>
int AccelModelTemplate<DeviceType, SolverType>::SetSomeParametersFromFile(
    int n, const std::string& filename, std::string& errorMessage)
{
    return device->SetSomeParametersFromFile(n, filename, errorMessage, dataFolder);
};
template <class DeviceType, class SolverType>
int AccelModelTemplate<DeviceType, SolverType>::SetSomeSolverFileName(const std::string& solverName,
                                                                      int fileNumber, int io,
                                                                      const std::string& filename,
                                                                      std::string& errorMessage)
{
    std::string filenameL = filename;
    //	if (!io)
    //	filenameL = dataFolder + filename;
    return solver->SetSomeSolverFileName(solverName, fileNumber, io, filenameL, errorMessage);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSomeSolverFileName(
    const std::string& solverName, std::vector<std::string>& input,
    std::vector<std::string>& output, std::vector<std::string>& filesInDescr,
    std::vector<std::string>& filesOutDescr)
{
    solver->GetSomeSolverFileName(solverName, input, output, filesInDescr, filesOutDescr);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSolverFileProps(const std::string& solverName,
                                                                    int                fileNumber,
                                                                    std::string&       ext,
                                                                    std::string&       prop)
{
    solver->GetSolverFileProps(solverName, fileNumber, ext, prop);
};
template <class DeviceType, class SolverType>
std::vector<std::string> AccelModelTemplate<DeviceType, SolverType>::GetSomeFileName()
{
    return device->GetSomeFileName();
};

template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetBrouseFlags(std::vector<std::string>& brouse,
                                                                std::vector<std::string>& names)
{
    return device->GetBrouseFlags(brouse, names);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::AddFlow()
{
    device->AddFlow();
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetParametersAccelFlow(
    std::vector<std::string>& keys, std::vector<double>& p, int flowNumber)
{
    device->GetLinacFlows()[flowNumber]->GetParametersAccelFlow(keys, p);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SetParametersAccelFlow(
    const std::vector<double>& params, int flowNumber)
{
    device->GetLinacFlows()[flowNumber]->SetParametersAccelFlow(params);
    device->GetLinacFlows()[flowNumber]->GenerateParticlesAccel();
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetMainAccelParameters(
    std::vector<std::string>& keysO, std::vector<double>& p)
{
    return device->GetMainAccelParameters(keysO, p);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetMainAccelParametersFlags(
    std::vector<std::string>& keysO, std::vector<double>& p)
{
    return device->GetMainAccelParametersFlags(keysO, p);
};

template <class DeviceType, class SolverType>
std::vector<double> AccelModelTemplate<DeviceType, SolverType>::GetMainAccelParameterCalculated()
{
    return device->GetMainAccelParameterCalculated();
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SetMainAccelParameters(
    const std::vector<double>& in)
{
    device->SetMainAccelParameters(in);
    device->TranslateParameters();
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SetSolverParameters(const std::string& solverName,
                                                                     const std::vector<double>& in)
{
    solver->SetSolverParameters(solverName, in);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSolverParameters(
    const std::string& solverName, std::vector<std::string>& keysO, std::vector<double>& p)
{
    solver->GetSolverParameters(solverName, keysO, p);
};

template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SetSolverParametersFlags(
    const std::string& solverName, const std::vector<double>& in)
{
    solver->SetSolverParametersFlags(solverName, in);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSolverParametersFlags(
    const std::string& solverName, std::vector<std::vector<std::string>>& keysO,
    std::vector<double>& p)
{
    solver->GetSolverParametersFlags(solverName, keysO, p);
};

template <class DeviceType, class SolverType>
std::vector<std::string> AccelModelTemplate<DeviceType, SolverType>::GetSolversNames()
{
    return solver->GetSolversNames();
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::GetSectionParameters(
    std::vector<std::vector<std::string>>& keysO, std::vector<std::vector<double>>& p)
{
    device->GetSectionParameters(keysO, p);
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SetAccelSectionParameters(
    int currentAccelSection, const std::vector<double>& params)
{
    device->SetAccelSectionParameters(currentAccelSection, params);
    device->CreateSequences();
};

template <class DeviceType, class SolverType>
AccelModelTemplate<DeviceType, SolverType>::AccelModelTemplate(const std::string& dataFolderIn)
{
    device     = std::shared_ptr<DeviceType>(new DeviceType());
    solver     = std::shared_ptr<SolverType>(new SolverType(dataFolderIn));
    dataFolder = dataFolderIn;
};

/*template <class DeviceType, class SolverType>
AccelModelTemplate<DeviceType, SolverType>::~AccelModelTemplate()
{
        //delete device;
        //delete solver;
}*/
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::SaveData(std::string fileName)
{
    std::ofstream                   ofs(fileName.c_str(), std::ios::out | std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << outputData;
};
template <class DeviceType, class SolverType>
void AccelModelTemplate<DeviceType, SolverType>::LoadData(std::string fileName)
{
    std::ifstream                   ifs(fileName.c_str(), std::ios::in | std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ia >> outputData;
};