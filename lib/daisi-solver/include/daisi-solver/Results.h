#ifndef RESULTS_H
#define RESULTS_H

#include "Geom.h"
#include "armadillo"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>

class SimulationData
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    std::vector<std::vector<std::string>>        dataFlags;
    std::vector<float>                           XData;
    std::vector<std::vector<std::vector<float>>> YDataFlow;
    std::vector<std::vector<std::vector<float>>> YDataElectrode;
    SimulationData();
    SimulationData(std::vector<int> flowsStyles, int nElectrodes, int type);
    template <class deviceType>
    void addDataPIC(const std::shared_ptr<deviceType>& device);

    void reset();
    void AddFlags(std::vector<std::string> dataFlagsIn);
    std::vector<std::vector<std::string>> GetFlags();

    template <class deviceType>
    void addData(deviceType& device, float error, int flag){
        /*		XData.push_back(device->GetFlow(0)->GetDynamicsData(0)->Time /
           (1e-9*commtools::LIGHT_VELOCITY()));

                        YData[0].push_back(error);

                        std::vector<float> emittance(3);


                        for (int i = 0; i < device->GetNumberParticlesFlows(); i++)
                        {
                                YDataFlow[i][0].push_back(device->GetFlow(i)->GetNumberOfParticles());
                                YDataFlow[i][1].push_back(0);
                                YDataFlow[i][2].push_back(0);

                                device->GetFlow(i)->GetRmsEmittances(emittance);

                                YDataFlow[i][3].push_back(emittance[0]);
                                YDataFlow[i][4].push_back(emittance[1]);
                                YDataFlow[i][5].push_back(emittance[2]);

                        }*/

        /*if (flag == 0)
        {
                for (int i = 0; i < device->GetNumberParticlesFlows(); i++)
                {
                        YDataFlow[i][3 * i +
1].push_back(device->GetFlow(i)->GetNumberOfParticles()); YDataFlow[i][3 * i +
2].push_back(device->GetFlow(i)->GetEmitterDevice()->GetEmissionCurrent()); float Eav =
device->GetFlow(i)->GetEmitterDevice()->getErAverage(); YDataFlow[i][3 * i + 3].push_back(Eav);
                }
        }
        if (flag == 1)
        {
                for (int i = 0; i < device->GetNumberParticlesFlows(); i++)
                {
                        YData[i][0].push_back(device->GetFlow(i)->GetNumberOfParticles());
                        YData[i][1].push_back(0);
                        YData[i][2].push_back(0);

                        device->GetFlow(i)->GetRmsEmittances(emittance);

                        YData[i][3].push_back(emittance[0]);
                        YData[i][4].push_back(emittance[1]);
                        YData[i][5].push_back(emittance[2]);

                }
        }
//		{ "N particles(t)", "I(t)", "Enorm_av(t)"};*/

    }
    void addData(int step, double bs);
};
class DynamicsData
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    std::vector<std::vector<unsigned int>> saveIndexes;
    std::vector<std::vector<int>>          writeIndexesPerThread;
    std::vector<std::vector<int>>          writeIndexes;
    std::string                            tag;

    int sizeElement;

  public:
    void addDataDyn(int index, float T, const std::vector<float>& X);
    void Init(int dataSize, float massIn, float chargeIn, int sizeElementIn, int NumberThreades,
              float lambdaIn);
    void InitL(int dataSize, float massIn, float chargeIn, int sizeElementIn, int NumberThreades,
               float lambdaIn, float EmIn);

    void AddSavedTraces(std::vector<unsigned int> saveIndexesIn, int threadNumber);
    void AddBlock(int saveInd1, int saveInd2, int threadNumber, int nBlocks, int blockNumber);

    std::vector<double>              TimeArrayAdd;
    std::vector<std::vector<double>> dataAdd;

    void SetAdditionalData(const std::vector<double>&        TimeArrayAddIn,
                           std::vector<std::vector<double>>& dataAddIn);
    void InitAdd(int dataSize, int nParticles);
    void SetDataAdd(std::vector<void*> dataIn, float Time);

    std::vector<std::vector<double>>              TimeArray;
    std::vector<std::vector<std::vector<double>>> data;

    std::vector<std::vector<double>>              TimeArray1;
    std::vector<std::vector<std::vector<double>>> data1;

    float            mass;
    float            charge;
    std::vector<int> steps;
    float            lambda;
    float            emittanceX;
    float            emittanceY;
    //	std::vector <ParticleTrace> data;

    DynamicsData();
    float StartTime();
    void SetEmptyPlaces(const std::vector<unsigned int>& EmptyPlaces, int threadNumber);
    void SetRemove(std::vector<unsigned int> EmptyPlaces, int threadNumber);
    void SetData(std::vector<void*> dataIn, float Time, int threadNumber, int saveParam, int flag);
    void Init(int dataSize, float massIn, float chargeIn, int sizeElementIn, int numberSaveTraces,
              int NumberThreades, float lambdaIn);
    void SetRemovePTI(std::vector<unsigned int> EmptyPlaces, int threadNumber);
    void AddBlock(std::vector<unsigned int> saveIndexesIn, int threadNumber, int nBlocks,
                  int blockNumber);
};

#ifdef SIM
class lineplot
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& flag;
        ar& line;
        ar& PlotTypeFlag;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& flag;
        ar& line;
        ar& PlotTypeFlag;
    }

  public:
    std::vector<std::string> flag;
    DGeo::Edge<double>       line;
    int                      PlotTypeFlag;
    lineplot(){

    };
    lineplot(std::vector<std::string> flagIn, double x1, double y1, double x2, double y2,
             int PlotTypeFlagIn)
    {
        flag          = flagIn;
        line.point2.y = y2;
        line.point2.x = x2;
        line.point2.z = 0;
        line.point1.y = y1;
        line.point1.x = x1;
        line.point1.z = 0;
        PlotTypeFlag  = PlotTypeFlagIn;
    };
    lineplot(std::vector<std::string> flagIn, double x1, double y1, double z1, double x2, double y2,
             double z2, int PlotTypeFlagIn)
    {
        flag          = flagIn;
        line.point2.y = y2;
        line.point2.x = x2;
        line.point2.z = z2;
        line.point1.y = y1;
        line.point1.x = x1;
        line.point1.z = z1;
        PlotTypeFlag  = PlotTypeFlagIn;
    };
};
#endif

class SynchrotronDevice;
class SimulationDataAccel : public DynamicsData
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    void SetDataNucl(arma::mat& dx, arma::mat& dy, arma::mat& dx0, arma::mat& dy0, arma::vec& dIx,
                     arma::vec& dIy, arma::vec& sx, arma::vec& sy);
    void addAccelElemetsDescription(std::shared_ptr<SynchrotronDevice>& device);
    std::vector<std::vector<float>>              XData;
    std::vector<std::vector<std::vector<float>>> YData;
    std::vector<std::vector<std::vector<float>>> YData1;

    std::vector<std::string>                                   dataFlags;
    std::string                                                tag;
    std::vector<std::vector<double>>                           props;
    std::vector<std::string>                                   names;
    std::vector<std::vector<std::vector<std::vector<double>>>> madXData;
    SimulationDataAccel(){}
    SimulationDataAccel(const std::vector<std::string>& dataFlagsInput, const std::string& tagIn,
                        int size, const std::vector<int>& Ysizes);
    SimulationDataAccel(const std::vector<std::string>& dataFlagsInput, const std::string& tagIn,
                        int size1, int size2);
    SimulationDataAccel(const std::vector<std::string>& dataFlagsInput, const std::string& tagIn)
        : dataFlags(dataFlagsInput), tag(tagIn) {}
    void addData(int position, float XDataIn, const std::vector<float>& YData);
    void addData(int position, float XDataIn, float YData);
    void initSynchrotron(int dataSize, int Ncircle, const std::vector<int>& ind,
                         const std::vector<int>& indP);
    void setDataSynchrotron(int i, int circle, const std::vector<int>& ind,
                            const std::vector<int>& indP, std::vector<void*> dataIn);
    void addData(int position, float XDataIn, const std::vector<double>& YDataIn);
};
#endif