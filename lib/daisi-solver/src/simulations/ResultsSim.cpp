#include "Results.h"
#include "DataTypes.h"
#include "Device.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "FlagStringsSolver.h"
#include "Particle.h"
#include "ParticlesFlow.h"
#include <Constants.h>

template void
SimulationData::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                      const unsigned int file_version);
template void
SimulationData::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                      const unsigned int file_version) const;

template <class Archive>
void SimulationData::save(Archive& ar, const unsigned int) const
{
    ar& XData;
    ar& YDataFlow;
    ar& YDataElectrode;
    ar& dataFlags;
}
template <class Archive>
void SimulationData::load(Archive& ar, const unsigned int)
{
    ar& XData;
    ar& YDataFlow;
    ar& YDataElectrode;
    ar& dataFlags;
}
void SimulationData::reset(){
#pragma ivdep
    /*for (int i = 0; i < YData.size(); i++)
            YData[i].clear();
    XData.clear();*/
};
void SimulationData::AddFlags(std::vector<std::string> dataFlagsIn){
    /*for (int i = 0; i < dataFlagsIn.size(); i++)
            dataFlags.push_back(dataFlagsIn[i]);

    YData.resize(dataFlags.size());*/
};
std::vector<std::vector<std::string>> SimulationData::GetFlags()
{
    return dataFlags;
};
void SimulationData::addData(int step, double bs)
{
    // XData.push_back(step);
    // YData[0].push_back(bs);
}

template void
DynamicsData::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                         const unsigned int file_version);
template void
DynamicsData::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                         const unsigned int file_version);

template <class Archive>
void DynamicsData::save(Archive& ar, const unsigned int) const
{
    ar& TimeArray;
    ar& data;
    ar& mass;
    ar& charge;
    ar& lambda;
    ar& TimeArrayAdd;
    ar& dataAdd;
    ar& tag;
}
template <class Archive>
void DynamicsData::load(Archive& ar, const unsigned int)
{
    ar& TimeArray;
    ar& data;
    ar& mass;
    ar& charge;
    ar& lambda;
    ar& TimeArrayAdd;
    ar& dataAdd;
    ar& tag;
}

void DynamicsData::SetAdditionalData(const std::vector<float>&        TimeArrayAddIn,
                                     std::vector<std::vector<float>>& dataAddIn)
{
    TimeArrayAdd = TimeArrayAddIn;
    dataAdd      = dataAddIn;
};
void DynamicsData::DynamicsData::InitAdd(int dataSize, int nParticles)
{
    TimeArray1.resize(nParticles);
    data1.resize(dataSize);

    for (int i = 0; i < dataSize; i++)
    {
        data1[i].resize(nParticles);
    }
}

template void
SimulationDataAccel::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void
SimulationDataAccel::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                                const unsigned int file_version);

template <class Archive>
void SimulationDataAccel::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<DynamicsData>(*this);
    ar& XData;
    ar& YData;
    ar& dataFlags;
    ar& tag;
}
template <class Archive>
void SimulationDataAccel::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<DynamicsData>(*this);
    ar& XData;
    ar& YData;
    ar& dataFlags;
    ar& tag;
}

SimulationData::SimulationData(){

};
SimulationData::SimulationData(std::vector<int> flowsStyles, int nElectrodes, int type)
{
    dataFlags.resize(flowsStyles.size() + 1);
    YDataFlow.resize(flowsStyles.size());

    YDataElectrode.resize(nElectrodes);
    for (int i = 0; i < nElectrodes; i++)
        YDataElectrode[i].resize(1);

    for (int i = 0; i < flowsStyles.size(); i++)
    {

        if (flowsStyles[i] < 5)
        {
            dataFlags[i] = flagStringsSolver::simulationDataNamesFlowPIC;
            YDataFlow[i].resize(flagStringsSolver::simulationDataNamesFlowPIC.size());
        }
        else
        {
            dataFlags[i] = flagStringsSolver::simulationDataNamesFlowPICAccel;
            YDataFlow[i].resize(flagStringsSolver::simulationDataNamesFlowPICAccel.size());
        }
    };
    dataFlags.back() = flagStringsSolver::simulationDataNamesElectrode;
};
template void
SimulationData::addDataPIC<device2daxsdouble>(const std::shared_ptr<device2daxsdouble>& device);
template void
SimulationData::addDataPIC<device2daxsfloat>(const std::shared_ptr<device2daxsfloat>& device);
template void
SimulationData::addDataPIC<device2ddouble>(const std::shared_ptr<device2ddouble>& device);
template void
SimulationData::addDataPIC<device2dfloat>(const std::shared_ptr<device2dfloat>& device);

template void
SimulationData::addDataPIC<device2dpolardouble>(const std::shared_ptr<device2dpolardouble>& device);
template void
SimulationData::addDataPIC<device2dpolarfloat>(const std::shared_ptr<device2dpolarfloat>& device);

template void
SimulationData::addDataPIC<device3dExtrfloat>(const std::shared_ptr<device3dExtrfloat>& device);
template void
SimulationData::addDataPIC<device3dExtrdouble>(const std::shared_ptr<device3dExtrdouble>& device);

template void
SimulationData::addDataPIC<device3ddouble>(const std::shared_ptr<device3ddouble>& device);

template <class deviceType>
void SimulationData::addDataPIC(const std::shared_ptr<deviceType>& device)
{
    XData.push_back(device->GetFlow(0)->GetDynamicsData(0)->Time /
                    (1e-9 * LIGHT_VELOCITY()));
    for (int i = 0; i < device->GetNumberParticlesFlows(); i++)
    {
        if (device->GetFlow(i)->GetDistributionStyle() < 5)
        {
            YDataFlow[i][0].push_back(float(device->GetFlow(i)->GetNumberOfParticles()));
            YDataFlow[i][1].push_back(device->GetFlow(i)->GetEmitterDevice()->GetEmissionCurrent());
            YDataFlow[i][2].push_back(device->GetFlow(i)->GetEmitterDevice()->getErAverage());
        }
        else
        {
            std::vector<float> emittance;
            device->GetFlow(i)->GetRmsEmittances(emittance);
            YDataFlow[i][0].push_back(float(device->GetFlow(i)->GetNumberOfParticles()));
            YDataFlow[i][1].push_back(emittance[0]);
            YDataFlow[i][2].push_back(emittance[1]);
            YDataFlow[i][3].push_back(emittance[2]);
        };
    };

    for (int i = 0; i < device->GetElectrodes().size(); i++)
        YDataElectrode[i][0].push_back(device->GetElectrodes()[i]->GetCurrent());
}