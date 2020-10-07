#include "Results.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"

#include <Constants.h>

SimulationDataAccel::SimulationDataAccel(const std::vector<std::string>& dataFlagsInput,
                                         const std::string& tagIn, int size1, int size2)
{
    tag       = tagIn;
    dataFlags = dataFlagsInput;
    XData.resize(size1);

    for (int i = 0; i < size1; i++)
    {
        XData[i].resize(size2);
    }
}

SimulationDataAccel::SimulationDataAccel(const std::vector<std::string>& dataFlagsInput,
                                         const std::string& tagIn, int size,
                                         const std::vector<int>& Ysizes)
{
    tag       = tagIn;
    dataFlags = dataFlagsInput;
    YData.resize(size);
    XData.resize(size);

    for (int i = 0; i < size; i++)
    {
        YData[i].resize(Ysizes[i]);
    }
}

void SimulationDataAccel::addData(int position, float XDataIn, const std::vector<float>& YDataIn)
{
    XData[position].push_back(XDataIn);
    YData[position].push_back(YDataIn);
}

void SimulationDataAccel::addData(int position, float XDataIn, const std::vector<double>& YDataIn)
{
    XData[position].push_back(XDataIn);
    YData[position].push_back({});
    for (int i = 0; i < YDataIn.size(); i++)
        YData[position].back().push_back(YDataIn[i]);
}

void SimulationDataAccel::addData(int position, float XDataIn, float YDataIn)
{
    XData[position].push_back(XDataIn);
    YData[position][0].push_back(YDataIn);
}

DynamicsData::DynamicsData()
{
    sizeElement = 0;
}

float DynamicsData::StartTime()
{
    if (TimeArray.size() == 0)
        return -1;

    if (TimeArray[0].size() == 0)
        return -1;

    return TimeArray[0][0];
}

void DynamicsData::SetEmptyPlaces(const std::vector<unsigned int>& EmptyPlaces, int threadNumber)
{
    int n = int(saveIndexes[threadNumber].size());
    for (int k = 0; k < EmptyPlaces.size(); k++)
    {
        for (int j = 0; j < n; j++)
        {
            if (saveIndexes[threadNumber][j] == EmptyPlaces[k])
            {
                saveIndexes[threadNumber].erase(saveIndexes[threadNumber].begin() + j);
                writeIndexes[threadNumber].erase(writeIndexes[threadNumber].begin() + j);
                if (n != 0)
                    n--;
                if (k != 0)
                    k--;
                if (j != 0)
                    j--;
            }
        }
    }
}

void DynamicsData::SetRemove(std::vector<unsigned int> EmptyPlaces, int threadNumber)
{
    int n  = int(saveIndexes[threadNumber].size());
    //int km = EmptyPlaces.size();

    n = int(saveIndexes[threadNumber].size());

    for (int j = 0; j < n; j++)
    {
        int s = 0;
        for (int k = 0; k < EmptyPlaces.size(); k++)
        {
            if (EmptyPlaces[k] < saveIndexes[threadNumber][j])
                s++;
        }
        saveIndexes[threadNumber][j] = saveIndexes[threadNumber][j] - s;
    }

    /*for (int k = 0; k < EmptyPlaces.size(); k++)
    {
            for (int j = 0; j < n; j++)
            {
                    if (EmptyPlaces[k] < saveIndexes[j])
                    {
                            for (int i = j; i < saveIndexes.size(); i++)
                                    saveIndexes[i]--;
                            break;
                    };
            }
    }*/
}

void DynamicsData::SetRemovePTI(std::vector<unsigned int> EmptyPlaces, int threadNumber)
{
    int n  = int(saveIndexes[threadNumber].size());
    int km = EmptyPlaces.size();
    for (int k = 0; k < km; k++)
    {
        for (int j = 0; j < n; j++)
        {
            if (saveIndexes[threadNumber][j] == EmptyPlaces[k])
            {
                writeIndexes[threadNumber][j] = -1;
                saveIndexes[threadNumber][j]  = -1;
            }
        }
    }

    if (saveIndexes[threadNumber].size() != 0)
    {
        n = int(writeIndexes[threadNumber].size());
        for (int k = 0; k < n; k++)
        {
            if (writeIndexes[threadNumber][k] == -1)
            {
                writeIndexes[threadNumber].erase(writeIndexes[threadNumber].begin() + k);
                saveIndexes[threadNumber].erase(saveIndexes[threadNumber].begin() + k);
                n--;
                k--;
            }
        }
    }
    std::vector<int> diff(saveIndexes[threadNumber].size());
    for (int i  = 0; i < saveIndexes[threadNumber].size(); i++)
        diff[i] = 0;

    for (int j = 0; j < EmptyPlaces.size(); j++)
    {
        for (int i = 0; i < saveIndexes[threadNumber].size(); i++)
        {
            if (saveIndexes[threadNumber][i] > EmptyPlaces[j])
                diff[i]++;
        }
    }

    for (int i                       = 0; i < saveIndexes[threadNumber].size(); i++)
        saveIndexes[threadNumber][i] = saveIndexes[threadNumber][i] - diff[i];
}

void DynamicsData::SetDataAdd(std::vector<void*> dataIn, float Time)
{
    if (sizeElement == 0)
        return;

    float tmp;

    for (int j = 0; j < dataIn.size(); j++)
    {
        for (int i = 0; i < data1[j].size(); i++)
        {
            if (sizeElement == 8)
            {
                tmp = *((double*)((char*)dataIn[j] + i * sizeElement));
                data1[j][i].push_back(tmp);
            }
        }
    }

    for (int i = 0; i < data1[0].size(); i++)
    {
        TimeArray1[i].push_back(Time / (1e-9));
    }
}

void DynamicsData::SetData(std::vector<void*> dataIn, float Time, int threadNumber, int saveParam,
                           int flag)
{
    if (sizeElement == 0)
        return;
    //int n = int(saveIndexes[threadNumber].size());

    float tmp;
    for (int j = 0; j < dataIn.size(); j++)
    {
        for (int i = 0; i < saveIndexes[threadNumber].size(); i++)
        {
            if (flag == 0 || 0 == steps[writeIndexes[threadNumber][i]] % saveParam)
            {
                if (sizeElement == 4)
                {
                    tmp =
                        *((float*)((char*)dataIn[j] + saveIndexes[threadNumber][i] * sizeElement));
                    data[j][writeIndexes[threadNumber][i]].push_back(tmp);
                }
                if (sizeElement == 8)
                {
                    tmp =
                        *((double*)((char*)dataIn[j] + saveIndexes[threadNumber][i] * sizeElement));

                    data[j][writeIndexes[threadNumber][i]].push_back(tmp);
                }
            }
        }
    }

    for (int i = 0; i < saveIndexes[threadNumber].size(); i++)
    {
        if (flag == 0 || 0 == steps[writeIndexes[threadNumber][i]] % saveParam)
        {
            TimeArray[writeIndexes[threadNumber][i]].push_back(Time / (1e-9));
        }
        if (flag == 1)
            steps[writeIndexes[threadNumber][i]]++;
    }
}

void DynamicsData::Init(int dataSize, float massIn, float chargeIn, int sizeElementIn,
                        int numberSaveTraces, int NumberThreades, float lambdaIn)
{
    if (numberSaveTraces <= 0)
        return;

    sizeElement = sizeElementIn;
    mass        = massIn;
    charge      = chargeIn;

    data.resize(dataSize);
    lambda        = lambdaIn;
    int perThread = numberSaveTraces / NumberThreades;

    writeIndexesPerThread.resize(NumberThreades);
    saveIndexes.resize(NumberThreades);

    for (int i = 0; i < NumberThreades; i++)
    {
        int i1 = perThread * i;
        int i2 = perThread * i + perThread;

        if (i == NumberThreades - 1)
            i2 = numberSaveTraces;

        for (int j = i1; j < i2; j++)
            writeIndexesPerThread[i].push_back(j);
    }

    TimeArray.resize(numberSaveTraces);
    for (int i = 0; i < dataSize; i++)
    {
        data[i].resize(numberSaveTraces);
    }
    writeIndexes = writeIndexesPerThread;
}

void DynamicsData::InitL(int dataSize, float massIn, float chargeIn, int sizeElementIn,
                         int NumberThreades, float lambdaIn, float EmIn)
{
    sizeElement = sizeElementIn;
    mass        = massIn;
    charge      = chargeIn;
    lambda      = lambdaIn;
    data.resize(dataSize);
    //	emittance = EmIn;

    writeIndexes.resize(NumberThreades);
    saveIndexes.resize(NumberThreades);
}

void DynamicsData::AddSavedTraces(std::vector<unsigned int> saveIndexesIn, int threadNumber)
{
    saveIndexes[threadNumber].insert(saveIndexes[threadNumber].end(), saveIndexesIn.begin(),
                                     saveIndexesIn.end());

    int currentTraces = TimeArray.size();
    int addTraces     = saveIndexesIn.size();
    TimeArray.resize(currentTraces + addTraces);
    steps.resize(currentTraces + addTraces);

    for (int i = 0; i < data.size(); i++)
    {
        data[i].resize(currentTraces + addTraces);
    }
    for (int i = 0; i < addTraces; i++)
    {
        steps[currentTraces + i] = 0;
        writeIndexes[threadNumber].push_back(i + currentTraces);
    }
}

void DynamicsData::Init(int dataSize, float massIn, float chargeIn, int sizeElementIn,
                        int NumberThreades, float lambdaIn)
{
    sizeElement = sizeElementIn;
    mass        = massIn;
    charge      = chargeIn;
    lambda      = lambdaIn;
    data.resize(dataSize);

    writeIndexes.resize(NumberThreades);
    saveIndexes.resize(NumberThreades);
}

void DynamicsData::AddBlock(std::vector<unsigned int> saveIndexesIn, int threadNumber, int nBlocks,
                            int blockNumber)
{
    int s = saveIndexesIn.size();

    int numberSaveTraces;
    numberSaveTraces = writeIndexesPerThread[threadNumber].size() / nBlocks;
    if (blockNumber == nBlocks - 1)
        numberSaveTraces =
            writeIndexesPerThread[threadNumber].size() - numberSaveTraces * (nBlocks - 1);

    if (numberSaveTraces < s)
    {
        for (int i = 0; i < s - numberSaveTraces; i++)
        {
            int k = rand() % saveIndexesIn.size();
            saveIndexesIn.erase(saveIndexesIn.begin() + k);
        }
    }
    steps.resize(numberSaveTraces);
    saveIndexes[threadNumber] = saveIndexesIn;
}

void DynamicsData::AddBlock(int saveInd1, int saveInd2, int threadNumber, int nBlocks,
                            int blockNumber)
{
    if (writeIndexesPerThread.size() == 0)
        return;

    std::vector<unsigned int> saveIndexesIn;

    for (int i = saveInd1; i < saveInd2; i++)
        saveIndexesIn.push_back(i);

    int s = saveIndexesIn.size();

    int numberSaveTraces;
    numberSaveTraces = writeIndexesPerThread[threadNumber].size() / nBlocks;
    if (blockNumber == nBlocks - 1)
        numberSaveTraces =
            writeIndexesPerThread[threadNumber].size() - numberSaveTraces * (nBlocks - 1);

    if (numberSaveTraces < s)
    {
        for (int i = 0; i < s - numberSaveTraces; i++)
        {
            int k = rand() % saveIndexesIn.size();
            saveIndexesIn.erase(saveIndexesIn.begin() + k);
        }
    }
    steps.resize(numberSaveTraces);
    saveIndexes[threadNumber] = saveIndexesIn;
}

void DynamicsData::addDataDyn(int index, float T, const std::vector<float>& X)
{
    if (index >= TimeArray.size())
    {
        TimeArray.push_back(std::vector<float>{});
        for (int i = 0; i < data.size(); i++)
            data[i].push_back(std::vector<float>{});
    };
    TimeArray[index].push_back(T);
    for (int i = 0; i < data.size(); i++)
        data[i][index].push_back(X[i]);
}

template <class Archive>
void SimulationDataAccel::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<DynamicsData>(*this);
    ar& XData;
    ar& YData;
    ar& dataFlags;
    ar& tag;
    ar& props;
    ar& names;
}
template <class Archive>
void SimulationDataAccel::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<DynamicsData>(*this);
    ar& XData;
    ar& YData;
    ar& dataFlags;
    ar& tag;
    ar& props;
    ar& names;
}

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

template void
DynamicsData::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                         const unsigned int file_version);
template void
DynamicsData::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                         const unsigned int file_version);

template void
SimulationDataAccel::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void
SimulationDataAccel::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                                const unsigned int file_version);

template void
SimulationDataAccel::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);

template void
SimulationDataAccel::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version) const;
template void
DynamicsData::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);

template void
DynamicsData::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version) const;
