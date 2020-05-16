#include "SynchrotronTools.h"
#include "../base/AccelFlow.h"
#include "Results.h"
#include "SynchrotronStrings.h"
#include "Tools.h"

#include <common_tools/constants.h>

char line[12500];
// char lineS[256];

void loadtrack(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData)
{
    outputData->Init(6, 1, 1, 8, 1, 1);
    std::vector<float>       tmp;
    std::vector<std::string> tmpStr;

    std::string        err;
    FILE*              fp    = fopen(fileName.c_str(), "r");
    bool               fhead = true;
    int                i     = 0;
    std::vector<float> X(6);
    float              S;
    int                index;
    int                turn;
    double             T, PT;
    float              currentAt = 0;
    bool               flagAt    = false;
    float              lastEnd   = 0;
    std::string        tmpS;

    while (fgets(line, 12500, fp))
    {
        if (strlen(line) && fhead)
        {
            switch (line[0])
            {
            case '@':
                break;
            case '$':
                fhead = false;
            };
        }
        else if (strlen(line) && !fhead)
        {
            if (line[0] != '#')
            {
                strsplit(line, " ", tmp, err);
                if (flagAt)
                {
                    lastEnd = tmp[8];
                    flagAt  = false;
                };
                outputData->addDataDyn(tmp[0] - 1, tmp[8] + currentAt,
                                       std::vector<float>(tmp.begin() + 2, tmp.begin() + 8));
            }
            else
            {
                tmpS = line;

                strsplit(tmpS, " ", tmpStr);
                if (tmpStr[5] == "end")
                {
                    currentAt = currentAt + lastEnd;
                    flagAt    = true;
                }
            };
        };
    }
    fclose(fp);
};
void loadtwiss(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData)
{
    std::vector<std::string> strs;
    FILE*                    fp    = fopen(fileName.c_str(), "r");
    bool                     fhead = true;
    int                      i     = 0;
    std::string              tmp;
    while (fgets(line, 12500, fp))
    {
        if (strlen(line) && fhead)
        {
            switch (line[0])
            {
            case '@':
                break;
            case '$':
                fhead = false;
            };
        }
        else if (strlen(line) && !fhead)
        {
            tmp = line;
            // strsplit(line, " \"", strs);
            strsplit(tmp, " \"", strs);

            for (int k = 0; k < 6; k++)
                outputData->addData(k, std::stod(strs[2]), std::stod(strs[3 + k]));
        };
    }
    fclose(fp);
};

void loadtwiss(const std::string& fileName, float& betaM, float& alphaM)
{
    std::vector<std::string> strs;
    FILE*                    fp    = fopen(fileName.c_str(), "r");
    bool                     fhead = true;
    int                      i     = 0;
    betaM                          = 0;
    alphaM                         = 0;
    while (fgets(line, 12500, fp))
    {
        if (strlen(line) && fhead)
        {
            switch (line[0])
            {
            case '@':
                break;
            case '$':
                fhead = false;
            };
        }
        else if (strlen(line) && !fhead)
        {
            strsplit(line, " \"", strs);

            if (std::abs(std::stod(strs[3])) > betaM)
                betaM = std::abs(std::stod(strs[3]));

            if (std::abs(std::stod(strs[6])) > betaM)
                betaM = std::abs(std::stod(strs[6]));

            if (std::abs(std::stod(strs[4])) > alphaM)
                alphaM = std::abs(std::stod(strs[4]));

            if (std::abs(std::stod(strs[7])) > alphaM)
                alphaM = std::abs(std::stod(strs[7]));
        };
    }
    fclose(fp);
}

void loadcm(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData)
{
    loadtrack(fileName, outputData);
    outputData->YData.resize(1);
    outputData->YData[0].push_back(outputData->TimeArray[0]);
    for (size_t k = 0; k < 4; k++)
    {
        outputData->YData[0].push_back(outputData->data[k][0]);
    }
    // outputData.back()->YData.resize(1);
    // outputData.back()->YData[0] = result;
}