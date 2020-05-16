#include <fstream>

#include <armadillo>

#include "../base/AccelFlow.h"
#ifdef WIN32
#include <windows.h>
#endif

#include <common_tools/constants.h>

#include "Results.h"
#include "SynchrotronDevice.h"
#include "SynchrotronSolver.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "Tools.h"

void SynchrotronSolver::SynchrotronSolversCorrEffStep(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    std::vector<std::string>& status, double sigma)
{

#ifdef WIN32
    TCHAR buffer[MAX_PATH];
    GetCurrentDirectory(MAX_PATH, &buffer[0]);
    std::wstring current = buffer;

    std::wstring workPath = current + L"/MADX";

    std::wstring configFile = workPath + L"/config.madx";

    std::wstring exePath = current + L"/MADX/madx-win64.exe " + configFile;

    STARTUPINFO         info = {sizeof(info)};
    PROCESS_INFORMATION processInfo;

    std::wstring fconF(inputFileNames[SynchrotronSolversCorrEffName][0].begin(),
                       inputFileNames[SynchrotronSolversCorrEffName][0].end());

    // int n = solverParameters[SynchrotronSolversCorrEffName]->find("Number of experiments");
    device->GetLinacFlows()[0]->SaveBetaToMadFile(std::string(current.begin(), current.end()) +
                                                  "/MADX/beta0.txt");

    /*if (!CopyFile(&fconF[0], &configFile[0], false))
    {
        errorMsg = errorMsg + "Unable to copy config.\n Please specify the config file path in the "
                              "\"solver parameters\"";
        return;
    }*/

    std::ofstream outfile;

    outfile.open(configFile, std::ofstream::out | std::ofstream::trunc);
    outfile.close();

    int nIter = solverParameters[SynchrotronSolversCorrEffName]->find("number of iterations");

    outfile.open(configFile, std::ios_base::app);
    outfile << "err_rel0_sigma = " << sigma << ";\n"
            << "err_rel1_sigma = " << 0 << ";\n"
            << "n_iter = " << nIter << ";\n\n\n";

    std::ifstream config;
    config.open(fconF, std::ios_base::in);

    outfile << config.rdbuf();

    outfile.close();
    config.close();

    if (CreateProcess(NULL, &exePath[0], NULL, NULL, TRUE, 0, NULL, &workPath[0], &info,
                      &processInfo))
    {

        WaitForSingleObject(processInfo.hProcess, INFINITE);
        CloseHandle(processInfo.hProcess);
        CloseHandle(processInfo.hThread);

        auto lambdaGetData = [&](const std::wstring file) {
            std::ifstream fS;
            fS.open(workPath + file, std::ifstream::in);

            int                ii = 0;
            std::vector<float> bx;
            std::vector<float> by;
            std::vector<float> xm;
            std::vector<float> ym;

            while (fS)
            {
                bx.push_back(0);
                by.push_back(0);
                xm.push_back(0);
                ym.push_back(0);

                fS >> bx[ii];
                fS >> by[ii];
                fS >> xm[ii];
                fS >> ym[ii];

                ii++;
            }
            bx.pop_back();
            by.pop_back();
            xm.pop_back();
            ym.pop_back();

            outputData.back()->XData.push_back(bx);
            outputData.back()->XData.push_back(by);
            outputData.back()->XData.push_back(xm);
            outputData.back()->XData.push_back(ym);

            fS.close();

            std::ofstream outfile;

            outfile.open(workPath + file, std::ofstream::out | std::ofstream::trunc);
            outfile.close();
        };
        lambdaGetData(L"/Tw_before_corr.txt");
        std::vector<float> xm = *(outputData.back()->XData.end() - 2);
        std::vector<float> ym = *(outputData.back()->XData.end() - 1);

        lambdaGetData(L"/Tw_corr.txt");

        std::vector<float> xmC = *(outputData.back()->XData.end() - 2);
        std::vector<float> ymC = *(outputData.back()->XData.end() - 1);

        for (size_t i = 0; i < xmC.size(); i++)
        {
            xmC[i] = xm[i] / xmC[i];
            ymC[i] = ym[i] / ymC[i];
        }

        outputData.back()->XData.push_back(xmC);
        outputData.back()->XData.push_back(ymC);

        return;
    }
#endif
    /*std::wstring resFile = workPath + L"/Ring-data.txt";
    for (int i = 0; i < n; i++)
    {
        if (CreateProcess(NULL, &exePath[0], NULL, NULL, TRUE, 0, NULL, &workPath[0], &info,
    &processInfo))
        {
            loadtwiss(std::string(resFile.begin(), resFile.end()), outputData.back()->XData[0][i],
                      outputData.back()->XData[1][i]);
        }
    }*/
}

void SynchrotronSolver::SynchrotronSolversCorrEff(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    std::vector<std::string>& status)
{
    double sig_start =
        solverParameters[SynchrotronSolversCorrEffName]->find("dip. comp. err. sigma start");

    double sig_end =
        solverParameters[SynchrotronSolversCorrEffName]->find("dip. comp. err. sigma end");

    int sig_n_steps = solverParameters[SynchrotronSolversCorrEffName]->find("number of steps");

    std::vector<double> sigs;
    if (sig_n_steps <= 1)
    {
        sigs.push_back(sig_start);
    }
    else
    {
        double d_s = (sig_end - sig_start) / (sig_n_steps - 1);
        for (int i = 0; i < sig_n_steps; i++)
        {
            sigs.push_back(sig_start + i * d_s);
        }
    }

    std::vector<std::string> Cur_SynchrotronEstFlags;
    for (int i = 0; i < sigs.size(); i++)
    {
        for (int j = 0; j < SynchrotronEstFlags.size(); j++)
        {
            Cur_SynchrotronEstFlags.push_back(SynchrotronEstFlags[j] + " " +
                                              std::to_string(sigs[i]));
        }
    }

    outputData.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(Cur_SynchrotronEstFlags, SynchrotronSolversCorrEffName)));

    for (int i = 0; i < sigs.size(); i++)
    {
        SynchrotronSolversCorrEffStep(device, progress, flagAbort, outputData, errorMsg, status,
                                      sig_start + sigs[i]);
    }
}