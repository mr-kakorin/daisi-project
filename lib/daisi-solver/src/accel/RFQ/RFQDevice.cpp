#include "RFQDevice.h"
#include "../base/AccelFlow.h"
#include "RFQStrings.h"
#include "Tools.h"

#include <common_tools/constants.h>

std::vector<std::vector<double>> params = {{107, -1.48, 1.9357, 1.12, 1.35709},
                                           {46, -1.4, 1, 1.13, 1},
                                           {143, -0.59, 1.5, 1.35, 5.57},
                                           {86, -0.55, 0.5, 1.4, 0.3}};
std::vector<double> RFQMainParametersDefault = {1.08e+8, 0.006, 4};

void RFQDevice::AddFlow()
{
    flows.push_back(std::shared_ptr<AccelFlow>(new AccelFlow(1)));
};

void RFQDevice::CreateSequences()
{
    Sequences.resize(3);
    double Phase0      = -commtools::PI() / 2;
    double Modulation0 = 1;
    double Phase1;
    double Modulation1;
    int    ncells;
    double degree;
    double degree1;
    Sequences[0].clear();
    Sequences[1].clear();
    Sequences[2].clear();

    int nCellMatcher = int(SectionParameters[0]->find("Number of cells"));

    for (int i = 0; i < nCellMatcher; i++)
    {
        Sequences[1].push_back(-commtools::PI() / 2);
        Sequences[0].push_back(1);
    };
    for (int section = 1; section < SectionParameters.size(); section++)
    {
        Phase1        = SectionParameters[section]->find("Output phase");
        Modulation1   = SectionParameters[section]->find("Output modulation");
        ncells        = int(SectionParameters[section]->find("Number of cells"));
        degree        = SectionParameters[section]->find("Appr. deg. phase");
        degree1       = SectionParameters[section]->find("Appr. deg. mod");
        double kMod   = (Modulation1 - Modulation0) / pow(double(ncells), degree1);
        double kPhase = (Phase1 - Phase0) / pow(double(ncells), degree);
        for (int i = 0; i < ncells; i++)
        {

            Sequences[0].push_back(Modulation0 + pow(double(i + 1), degree1) * kMod);
            Sequences[1].push_back(Phase0 + pow(double(i + 1), degree) * kPhase);
        };
        Modulation0 = Modulation1;
        Phase0      = Phase1;
    }
    Sequences[1].push_back(Sequences[1].back());

    double mD = SectionParameters[0]->find("Appr. deg.");

    double Rmax = SectionParameters[0]->find("Maximal Radius");

    double kMatcher =
        (Rmax - mainAccelParameters->find("Channel Radius")) / (pow(nCellMatcher, mD));

    Sequences[2].resize(nCellMatcher);

    for (int i = 0; i < nCellMatcher; i++)
    {
        Sequences[2][i] =
            mainAccelParameters->find("Channel Radius") + pow(nCellMatcher - i, mD) * kMatcher;
    }
};

bool RFQDevice::checkSequence()
{
    for (size_t i = 1; i < Sequences[0].size(); i++)
    {
        if (Sequences[0][i] < Sequences[0][i - 1])
        {
            return false;
        }
    }
    for (size_t i = 1; i < Sequences[1].size(); i++)
    {
        if (Sequences[1][i] < Sequences[1][i - 1])
        {
            return false;
        }
    }
    return true;
}

int RFQDevice::SetSomeParametersFromFile(int n, const std::string& filename,
                                         std::string& errorMessage, const std::string& folder)
{
    return 0;
}

void RFQDevice::TranslateParameters()
{
}

RFQDevice::RFQDevice()
{
    mainAccelParameters = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    for (int i = 0; i < RFQMainParameters.size(); i++)
        mainAccelParameters->insert(RFQMainParameters[i], RFQMainParametersDefault[i]);

    SectionParameters.resize(RFQMainParametersDefault[2] + 1);
    SectionParameters[0] = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    SectionParameters[0]->insert(RFQSectionNames0[0], 8);
    SectionParameters[0]->insert(RFQSectionNames0[1], 0.045);

    mainAccelParametersFlags = std::shared_ptr<myunsorted_map>(new myunsorted_map());

    for (int j = 1; j < RFQMainParametersDefault[2] + 1; j++)
    {
        SectionParameters[j] = std::shared_ptr<myunsorted_map>(new myunsorted_map());
        for (int i = 0; i < RFQSectionNames.size(); i++)
            SectionParameters[j]->insert(RFQSectionNames[i], params[j - 1][i]);
    }

    namesOfSequences = {"Modulations", "Sync phases", "Matcher radii"};
    CreateSequences();

    calcParametersNames = {"Cells lengths", "Minimal radii", "Minimal regular", "Acceleration Eff",
                           "Energy"};
    calcParameters.resize(calcParametersNames.size());

    calcMainParameters = std::shared_ptr<myunsorted_map>(new myunsorted_map());
};
