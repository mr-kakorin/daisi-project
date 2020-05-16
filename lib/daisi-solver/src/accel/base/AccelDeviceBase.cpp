#include "AccelDeviceBase.h"
#include "../base/AccelFlow.h"
#include "Tools.h"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

template class AccelDeviceBase<AccelFlow>;
template class AccelDeviceBase<SynchrotronFlow>;

template void AccelDeviceBase<AccelFlow>::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void AccelDeviceBase<AccelFlow>::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version);
template void AccelDeviceBase<SynchrotronFlow>::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void AccelDeviceBase<SynchrotronFlow>::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version);
/*AccelDeviceBase::~AccelDeviceBase()
{
        delete mainAccelParameters;
};*/
template <class FlowType>
std::vector<std::string> AccelDeviceBase<FlowType>::GetcalcMainParameters()
{
    return calcMainParameters->GetKeys();
};

template <class FlowType>
std::vector<double> AccelDeviceBase<FlowType>::GetMainAccelParameterCalculated()
{
    return calcMainParameters->GetValues();
};

template <class FlowType>
void AccelDeviceBase<FlowType>::checkInputParameters(std::string& errorMessage)
{
    for (int i = 0; i < mainAccelParameters->GetValues().size(); i++)
    {
        if (!mainAccelParameters->GetValues()[i])
        {
            errorMessage = errorMessage + "Zero main accelerator parameter.";
            return;
        };
    };
};
template <class FlowType> std::vector<double>& AccelDeviceBase<FlowType>::GetMainAccelParameters()
{
    return mainAccelParameters->GetValues();
};

template <class FlowType>
std::vector<std::vector<double>>& AccelDeviceBase<FlowType>::GetSequencesOfParameters()
{
    return Sequences;
};
template <class FlowType>
std::vector<std::vector<double>>& AccelDeviceBase<FlowType>::GetSequencesOfParametersCalc()
{
    return calcParameters;
};

template <class FlowType>
std::vector<std::string>& AccelDeviceBase<FlowType>::GetcalcParametersNames()
{
    return calcParametersNames;
};

template <class FlowType>
void AccelDeviceBase<FlowType>::GetSequencesOfParameters(int seq, std::vector<double>& sequence)
{
    sequence = Sequences[seq];
};

template <class FlowType>
void AccelDeviceBase<FlowType>::GetSequencesOfParametersCalc(int seq, std::vector<double>& sequence)
{
    sequence = calcParameters[seq];
};
template <class FlowType> std::vector<std::string>& AccelDeviceBase<FlowType>::GetnamesOfSequences()
{
    return namesOfSequences;
};

template <class FlowType>
void AccelDeviceBase<FlowType>::SetAccelSectionParameters(int currentAccelSection,
                                                          const std::vector<double>& params)
{
    SectionParameters[currentAccelSection]->SetValues(params);
};

template <class FlowType>
void AccelDeviceBase<FlowType>::GetSectionParameters(std::vector<std::vector<std::string>>& keysO,
                                                     std::vector<std::vector<double>>&      p)
{
    keysO.clear();
    p.clear();
    for (int i = 0; i < SectionParameters.size(); i++)
    {
        keysO.push_back(SectionParameters[i]->GetKeys());
        p.push_back(SectionParameters[i]->GetValues());
    }
};
template <class FlowType> std::vector<std::string>& AccelDeviceBase<FlowType>::GetSomeFileName()
{
    return files;
};

template <class FlowType>
std::vector<std::shared_ptr<FlowType>>& AccelDeviceBase<FlowType>::GetLinacFlows()
{
    return flows;
};
template <class FlowType>
void AccelDeviceBase<FlowType>::GetMainAccelParameters(std::vector<std::string>& keysO,
                                                       std::vector<double>&      p)
{
    keysO = mainAccelParameters->GetKeys();
    p     = mainAccelParameters->GetValues();
};
template <class FlowType>
void AccelDeviceBase<FlowType>::GetMainAccelParametersFlags(std::vector<std::string>& keysO,
                                                            std::vector<double>&      p)
{
    keysO = mainAccelParametersFlags->GetKeys();
    p     = mainAccelParametersFlags->GetValues();
};
template <class FlowType>
void AccelDeviceBase<FlowType>::SetMainAccelParameters(const std::vector<double>& in)
{
    mainAccelParameters->SetValues(in);
};
template <class FlowType> double AccelDeviceBase<FlowType>::GetParameter(const std::string& key)
{
    return mainAccelParameters->find(key);
};
template <class FlowType>
void AccelDeviceBase<FlowType>::GetBrouseFlags(std::vector<std::string>& brouse,
                                               std::vector<std::string>& names)
{
    brouse = filesExtensions;
    names  = filesNames;
};

template <class FlowType>
template <class Archive>
void AccelDeviceBase<FlowType>::save(Archive& ar, const unsigned int) const
{
    ar& mainAccelParameters;
    ar& SectionParameters;
    ar& flows;
    ar& filesNames;
    ar& filesExtensions;
    ar& files;
    ar& namesOfSequences;
    ar& Sequences;
    ar& calcParametersNames;
    ar& calcParameters;
    ar& calcMainParameters;
    ar& mainAccelParametersFlags;
};
template <class FlowType>
template <class Archive>
void AccelDeviceBase<FlowType>::load(Archive& ar, const unsigned int)
{
    ar& mainAccelParameters;
    ar& SectionParameters;
    ar& flows;
    ar& filesNames;
    ar& filesExtensions;
    ar& files;
    ar& namesOfSequences;
    ar& Sequences;
    ar& calcParametersNames;
    ar& calcParameters;
    ar& calcMainParameters;
    ar& mainAccelParametersFlags;
};
