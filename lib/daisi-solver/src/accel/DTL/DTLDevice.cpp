#include "DTLDevice.h"
#include "../base/AccelFlow.h"
#include "DTLStrings.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Tools.h"
#include <common_tools/constants.h>
std::vector<std::string> DTLSectionNames0 = {"Start phase", "Number of cells", "Max phase", "Appr. deg. phase (p)"};
std::vector<std::string> DTLSectionNames  = {"Number of cells", "Max phase", "Appr. deg. phase (p)"};
std::vector<std::vector<double>> params   = {{-1.5, 20, -1.57, 2}, {20, 1.57, 2}, {20, -1.57, 2}, {20, 1.57, 2},
                                           {20, -1.57, 2},       {20, 1.57, 2}, {20, -1.57, 2}};
// std::vector<double> L0 = { 0.02,0.02,0.02,0.02 };
std::vector<double>      DTLMainParametersDefault = {4.33e+8, 0.006, 7};
std::vector<std::string> DTLMainParameters        = {"Frequency, Hz", "Channel Radius", "Number of sections"};

void DTLDevice::AddFlow()
{
    flows.push_back(std::shared_ptr<AccelFlow>(new AccelFlow(2)));
};

void DTLDevice::CreateSequences()
{
    Sequences.resize(2);
    Sequences[0].clear();
    Sequences[1].clear();
    int SequenceMode = int(SectionParameters[0]->find("Number of cells"));
    switch (SequenceMode)
    {
    case -1:
    {
        std::string path1 = "C:/Users/Irina/Documents/Visual Studio "
                            "2015/Projects/DaisiShareNewMerge/DaisiClientFull/x64/Debug/Seq60.txt";
        std::ifstream ifs1;
        ifs1.open(path1, std::ios_base::in);
        double ph;
        while (ifs1.good())
        {
            ifs1 >> ph;
            Sequences[0].push_back(ph * 2 * commtools::PI() / 360);
        }
        break;
    }
    case -2:
    {
        int nCells = int(SectionParameters[1]->find("Number of cells"));
        for (int j = 0; j < nCells; j++)
            Sequences[0].push_back(SectionParameters[0]->find("Start phase"));
        Sequences[1].push_back(Sequences[0].size());

        break;
    }
    default:
    {
        std::vector<double> MaxPhase = {-commtools::PI() / 2, commtools::PI() / 2, -commtools::PI() / 2, commtools::PI() / 2,
                                        -commtools::PI() / 2, commtools::PI() / 2, -commtools::PI() / 2};
        int    ncells;
        double degree;
        double stphase;
        Sequences[0].clear();
        std::vector<double> Temp;
        std::vector<double> Temp1;
        // std::vector<double>coeff(2);
        // std::vector<double>L = L0;
        // double PointsPerCell = 10;

        ncells      = int(SectionParameters[1]->find("Number of cells"));
        degree      = SectionParameters[0]->find("Appr. deg. phase (p)");
        stphase     = SectionParameters[0]->find("Start phase"); // not used really
        int counter = 0;
        for (int k = ncells / 2; k < ncells + 1; k++)
        {
            Temp.push_back((-1) * MaxPhase[0] * pow(2 * (double(k) / double(ncells) - 0.5), degree) + MaxPhase[0]);
            Temp1.push_back((-1) * MaxPhase[0] * pow(2 * (double(k) / double(ncells) - 0.5), degree) + MaxPhase[0]);
        }
        /*
        for (int k = 0; k < ncells / 2; k++)
        {
        Temp1.pop_back();
        if (Temp1.back() < (stphase+0.1))
        {
        Sequences[0].push_back(Temp1.back());
        counter++;
        }
        }
        */
        // SectionParameters[0]->SetValue("Number of cells", ncells / 2);// +counter);
        for (int k = 1; k < Temp.size(); k++)
        {
            Sequences[0].push_back(Temp[k]);
        }
        // arma::mat A;//A*coeff=B
        // arma::vec B;
        // arma::vec x;
        // A << pow(5, degree) << pow(5, degree/2) << arma::endr << pow(ncells, degree) << pow(ncells, degree/2) <<
        // arma::endr<<-10<<1<< arma::endr;  B << MaxPhase[0]-stphase << (-1)*stphase<<0;  x = arma::solve(A, B);
        // std::copy(x.begin(), x.end(), coeff.begin());

        for (int section = 1; section < SectionParameters.size(); section++)
        {
            ncells = int(SectionParameters[section]->find("Number of cells"));
            degree = SectionParameters[section]->find("Appr. deg. phase (p)");

            Temp.clear();
            Temp1.clear();
            // A <<pow(ncells /2, degree)<< pow(ncells / 2, degree/2)<<arma::endr <<  pow(ncells, degree)<< pow(ncells,
            // degree/2) <<arma::endr;  B << MaxPhase[section]<<0;  x = arma::solve(A, B);  std::copy(x.begin(),
            // x.end(), coeff.begin());

            for (int k = ncells / 2; k < ncells + 1; k++)
            {
                Temp.push_back((-1) * MaxPhase[section] * pow(2 * (double(k) / double(ncells) - 0.5), degree) +
                               MaxPhase[section]);
                Temp1.push_back((-1) * MaxPhase[section] * pow(2 * (double(k) / double(ncells) - 0.5), degree) +
                                MaxPhase[section]);
            }
            int a = Temp.size();
            int b = Temp1.size();
            int c = Sequences[0].size();
            for (int k = 0; k < ncells / 2; k++)
            {
                Temp1.pop_back();
                Sequences[0].push_back(Temp1.back());
            }
            a = Temp.size();
            b = Temp1.size();
            c = Sequences[0].size();
            for (int k = 1; k < Temp.size(); k++)
            {
                Sequences[0].push_back(Temp[k]);
            }
            c = Sequences[0].size();
        }

        Sequences[1].push_back(Sequences[0].size());
        break;
    }
    }
};

int DTLDevice::SetSomeParametersFromFile(int n, const std::string& filename, std::string& errorMessage, const std::string& folder)
{
    return 0;
};

void DTLDevice::TranslateParameters(){

};

DTLDevice::DTLDevice()
{
    mainAccelParameters = std::shared_ptr<myunsorted_map>(new myunsorted_map());

    mainAccelParametersFlags = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    for (int i = 0; i < DTLMainParameters.size(); i++)
        mainAccelParameters->insert(DTLMainParameters[i], DTLMainParametersDefault[i]);

    if (SectionParameters.empty())
    {
        SectionParameters.resize(DTLMainParametersDefault[2]);

        for (int j = 0; j < DTLMainParametersDefault[2]; j++)
        {
            SectionParameters[j] = std::shared_ptr<myunsorted_map>(new myunsorted_map());
            if (j == 0)
            {
                for (int i = 0; i < DTLSectionNames0.size(); i++)
                    SectionParameters[j]->insert(DTLSectionNames0[i], params[j][i]);
            }
            else
            {
                for (int i = 0; i < DTLSectionNames.size(); i++)
                    SectionParameters[j]->insert(DTLSectionNames[i], params[j][i]);
            }
        }
    }

    namesOfSequences = {"Sync phases"};
    CreateSequences();

    calcParametersNames = {"Cells lengths", "Gap lengths", "Inner radii", "Outer radii", "Energy"};
    calcParameters.resize(calcParametersNames.size());

    calcMainParameters = std::shared_ptr<myunsorted_map>(new myunsorted_map());
};
