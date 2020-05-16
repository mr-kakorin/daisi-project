#define NOMINMAX

#include <fstream>

#include <common_tools/constants.h>
#include <notk/controller.hpp>

#include "../base/AccelFlow.h"
#include "Results.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "SynchrotronToolsSim.h"

#include "SynchrotronDevice.h"
#include "SynchrotronSolver.h"
#include "Tools.h"

#ifdef WIN32
#include <Windows.h>
#endif

double calcx(double x12, double x23, double s, double s1, double s2, double s3)
{
    double p1 = (s3 + s2 - 2 * s) / (s3 - s1);
    double p2 = (2 * s - s1 - s2) / (s3 - s1);
    return p1 * x12 + p2 * x23;
}

int findMaxel(const std::vector<double>& X, double x)
{
    std::vector<int> ind;
    for (int i = 0; i < X.size(); i++)
    {
        if (X[i] <= x)
            ind.push_back(i);
    };
    if (!ind.size())
        return -1;

    int    maxI = 0;
    double max  = ind[0];
    for (int i = 1; i < ind.size(); i++)
    {
        if (ind[i] > ind[max])
        {
            maxI = i;
            max  = ind[i];
        };
    }
    return maxI;
};

void tabl_ind3(std::vector<int>& ind, const std::vector<double>& X, double x)
{
    int i = findMaxel(X, x);
    int n = X.size();
    if (i != -1 && i < n - 1)
    {
        if (i > 0)
        {
            if (x <= X[i] + (X[i + 1] - X[i]) / 2)
                ind = {i - 1, i, i + 1};
        }
        if (i <= n - 3)
        {
            if (x >= X[i] + (X[i + 1] - X[i]) / 2)
                ind = {i, i + 1, i + 2};
        }
    };
};

void eqsolve(double& a, double& phi, double x1, double x2, double s1, double s2)
{
    arma::mat A(2, 2);

    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 1;
    A(1, 1) = -1;

    arma::vec b(2);
    b(0) = 2 * atan((x1 + x2) / (x1 - x2) * tan((s1 - s2) / 2));
    b(1) = s1 - s2;

    arma::vec y = arma::inv(A) * b;
    phi         = y(0) - s1;
    a           = x1 / std::sin(s1 + phi);
}

void SynchrotronSolver::ClosedOrbitCalculation(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{

    progress = 0;

    auto pt = readJSONFile(inputFileNames[SynchrotronSolversNameClosedOrbit][0]);

    if (!pt)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto errors = read_errors(*pt);

    if (!errors)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto seq_err = std::make_shared<OpticElementsSequence>(*device->GetOpticElementsSequence());

    seq_err->InsertErrors(*errors);

    seq_err->mtrans();

    auto result = calc_closed_orbit(seq_err, inputFileNames[SynchrotronSolversNameClosedOrbit][1],
                                    device->GetParameter("X aperture, m"),
                                    device->GetParameter("Y aperture, m"));

    if (!result.second.empty())
    {
        errorMsg = result.second;
        return;
    }

    device->GetLinacFlows()[0]->SetMassCenterVector(result.first.first, false);
    device->GetLinacFlows()[0]->SetMassCenterVector(result.first.second, true);

    progress = 1.0;
}

void SynchrotronSolver::DynamicsSimulationTwiss(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    bool flagSave)
{
    if (flagSave)
        outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
            SynchrotronTwissDynFlags, SynchrotronSolversNameTwiss, 6, {1, 1, 1, 1, 1, 1})));

    std::vector<arma::mat> xtwiss;
    std::vector<arma::mat> ytwiss;

    std::vector<float> mu_y;
    std::vector<float> mu_x;

    auto result = DynamicsSimulationTwissF<float>(device->GetLinacFlows()[0]->GetTwissVector(),
                                                  device->GetOpticElementsSequence(), 1);

    std::vector<float> S = result[0];

    for (int i = 0; i < device->GetOpticElementsSequence()->length(); i++)
    {
        if (!flagSave)
            continue;

        for (int k = 0; k < 2; k++)
            outputData.back()->addData(k, S[i], result[k + 1][i]);

        outputData.back()->addData(2, S[i], result[7][i]);

        for (int k = 0; k < 2; k++)
            outputData.back()->addData(k + 3, S[i], result[k + 4][i]);

        outputData.back()->addData(5, S[i], result[8][i]);
    };
    if (flagSave)
    {
        //   outputData.back()->addAccelElemetsDescription(device);
        /* savetwiss(outpuFilePaths[SynchrotronSolversNameTwiss][0], xtwiss, ytwiss,
                   device->GetLinacFlows()[0]->getParticleType(),
                   device->GetLinacFlows()[0]->getRestMassInGeV(),
                   device->GetOpticElementsSequence(), mu_x, mu_y);*/
    }
    progress = 1;
};

void SynchrotronSolver::DynamicsSimulationTwissMadx(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{
// std::wstring file(inputFileNames[SynchrotronSolversMADX][0].begin(),
// inputFileNames[SynchrotronSolversMADX][0].end());

#ifdef WIN32

    TCHAR buffer[MAX_PATH];
    GetCurrentDirectory(MAX_PATH, &buffer[0]);
    std::wstring current = buffer;

    std::wstring path2 = current + L"/MADX";

    std::wstring file = path2 + L"/config.madx";
    device->SaveMADXConfigFile(std::string(file.begin(), file.end()));

    std::wstring path1 = current + L"/MADX/madx-win64.exe " + file;

    STARTUPINFO         info = {sizeof(info)};
    PROCESS_INFORMATION processInfo;
    // UINT error = WinExec(RunKey.c_str(), 1);

    std::wstring fOpt(device->GetSomeFileName()[0].begin(), device->GetSomeFileName()[0].end());
    std::wstring fOptdest = path2 + L"/opt.opt";
    // fOptdest = current + L"/opt.opt";
    if (!CopyFile(&fOpt[0], &fOptdest[0], false))
    {
        errorMsg = errorMsg + "Unable to copy *.opt file.\n Please specify the *.opt file path in "
                              "the \"Accelerator parameters section\"";
        return;
    }

    std::wstring fLine(device->GetSomeFileName()[1].begin(), device->GetSomeFileName()[1].end());
    std::wstring flinedest = path2 + L"/line.line";
    // flinedest = current + L"/line.line";
    if (!CopyFile(&fLine[0], &flinedest[0], false))
    {
        errorMsg = errorMsg + "Unable to copy *.line file.\n Please specify the *.line file path "
                              "in the \"Accelerator "
                              "parameters section\"";
        return;
    }
    device->GetLinacFlows()[0]->SaveToMadFile(std::string(current.begin(), current.end()) +
                                              "/MADX/init.txt");
    device->GetLinacFlows()[0]->SaveBetaToMadFile(std::string(current.begin(), current.end()) +
                                                  "/MADX/beta0.txt");

    device->GetOpticElementsSequence()->SaveMADObsCommands(
        std::string(current.begin(), current.end()) + "/MADX/obs.txt");

    std::string align = std::string(current.begin(), current.end()) + "/MADX/align.txt";

    if (solverParametersFlags[SynchrotronSolversMADX][0] == 0)
    {
        device->GetOpticElementsSequence()->SaveMADAlignmentCommands(align);
    }
    else
    {
        FILE* fid;
        fopen_s(&fid, align.c_str(), "w");
        fclose(fid);
    }

    if (CreateProcess(NULL, &path1[0], NULL, NULL, TRUE, 0, NULL, &path2[0], &info, &processInfo))
    {
        WaitForSingleObject(processInfo.hProcess, INFINITE);
        CloseHandle(processInfo.hProcess);
        CloseHandle(processInfo.hThread);

        std::wstring f1 = path2 + L"/twiss_mad.txt";
        std::wstring f1dest(outpuFilePaths[SynchrotronSolversMADX][0].begin(),
                            outpuFilePaths[SynchrotronSolversMADX][0].end());

        if (!CopyFile(&f1[0], &f1dest[0], false))
            errorMsg = errorMsg + "Unable to copy twiss_mad.txt file";
        else
        {
            outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
                SynchrotronTwissDynFlags, SynchrotronSolversNameTwiss + " MADX", 8,
                {1, 1, 1, 1, 1, 1, 1, 1})));
            loadtwiss(outpuFilePaths[SynchrotronSolversMADX][0], outputData.back());
        };

        std::wstring f2 = path2 + L"/beam_one";
        std::wstring f2dest(outpuFilePaths[SynchrotronSolversMADX][1].begin(),
                            outpuFilePaths[SynchrotronSolversMADX][1].end());

        if (!CopyFile(&f2[0], &f2dest[0], false))
            errorMsg = errorMsg + "Unable to copy trackone file";
        else
        {
            outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
                SynchrotronDynFlags, SynchrotronSolversBeam + " MADX", 6, {1, 1, 1, 1, 1, 1})));
            loadtrack(outpuFilePaths[SynchrotronSolversMADX][1], outputData.back());

            std::wstring f3 = path2 + L"/cm_one";

            outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
                SynchrotronDynCMFlags, SynchrotronSolversNameCenterOfMass + " MADX")));

            loadcm(std::string(f3.begin(), f3.end()), outputData.back());
        };
    }
    else
    {
        errorMsg = "Unable to start process madx-win64.exe";
    };
        //	WaitForSingleObject(pi.hProcess, INFINITE);

#endif

    progress = 1;
};

void SynchrotronSolver::DynamicsSimulationBeam(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{

    if (!device->GetLinacFlows()[0]->GetnParticles())
    {
        errorMsg = "Number of particles should be > 0";
        return;
    };

    auto pt = readJSONFile(inputFileNames[SynchrotronSolversBeam][0]);

    if (!pt)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto errors = read_errors(*pt);

    if (!errors)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto seq_err = std::make_shared<OpticElementsSequence>(*device->GetOpticElementsSequence());

    seq_err->InsertErrors(*errors);

    seq_err->mtrans();

    outputData.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(SynchrotronDynFlags, SynchrotronSolversBeam)));

    int nParticles = device->GetLinacFlows()[0]->GetnParticles();
    device->GetLinacFlows()[0]->GetDynamicsAccel().Init();

    int ns = std::min(nParticles, int(solverParameters[SynchrotronSolversBeam]->find(
                                      "Number of visualized traces")));
    outputData.back()->Init(6, device->GetLinacFlows()[0]->GetMass(),
                            device->GetLinacFlows()[0]->GetCharge(), sizeof(double), ns, 1, 0);
    outputData.back()->AddBlock(0, ns, 0, 1, 0);
    int ncircle = device->GetParameter("Number of cicles");

    //������ ���������, �� ������� ����� � ����
    std::vector<int> ind = {0,  1,  3,   9,
                            19, 99, 299, int(device->GetOpticElementsSequence()->length() - 1)};
    //������ ������, ������� ����� � ����
    std::vector<int> indP = {0, 1, 2};

    outputData.back()->initSynchrotron(6, ncircle, ind, indP);

    double                          atC = 0;
    std::vector<std::vector<float>> OutParams(3);
    std::vector<float>              TimeOut;

    simulationBeam<float>(
        outputData.back()->TimeArray, outputData.back()->data, device->GetLinacFlows()[0], seq_err,
        device->GetParameter("Number of cicles"), progress, flagAbort, "ALL", 0, ns);

    OutParams[0].resize(outputData.back()->TimeArray[0].size());
    OutParams[1].resize(outputData.back()->TimeArray[0].size());

    std::fill(OutParams[0].begin(), OutParams[0].end(), device->GetParameter("X aperture, m"));
    std::fill(OutParams[1].begin(), OutParams[1].end(), device->GetParameter("Y aperture, m"));

    OutParams[2] =
        calculateTransmission(outputData.back()->data, device->GetParameter("X aperture, m"),
                              device->GetParameter("Y aperture, m"));

    outputData.back()->addAccelElemetsDescription(device);
    outputData.back()->SetAdditionalData(outputData.back()->TimeArray[0], OutParams);

    progress = 1;
    /*  for (int j = 0; j < ncircle; j++)
{
    double at = 0;
    for (int i = 0; i < device->GetOpticElementsSequence()->length(); i++)
    {
        updateBeamPositions(device->GetOpticElementsSequence()->Tx[i],
                            device->GetOpticElementsSequence()->Ty[i],
                            device->GetLinacFlows()[0]->GetDynamicsAccel());
        at = atC + device->GetOpticElementsSequence()->getParameter(i, "at") +
             device->GetOpticElementsSequence()->getParameter(i, "L");
        outputData.back()->SetData(device->GetLinacFlows()[0]->GetDynamicsAccel().GetData(),
                                   at * 1e-9, 0, 1, 0);
        outputData.back()->setDataSynchrotron(
            i, j, ind, indP, device->GetLinacFlows()[0]->GetDynamicsAccel().GetData());
        progress = double(j * device->GetOpticElementsSequence()->length() + i) /
                   double(ncircle * device->GetOpticElementsSequence()->length());
        OutParams[0].push_back(device->GetParameter("X aperture, m"));
        OutParams[1].push_back(device->GetParameter("Y aperture, m"));
        TimeOut.push_back(at);

        OutParams[2].push_back(
            device->GetLinacFlows()[0]->GetDynamicsAccel().CalculateTransmission(
                device->GetParameter("X aperture, m"), device->GetParameter("Y aperture,
m")));
    }
    atC = at;
}
//	outputData.back()->setDataSynchrotron(device->GetOpticElementsSequence()->length() -
1,
// ncircle-1,  ind, indP,
// device->GetLinacFlows()[0]->GetDynamicsAccel().GetData());
//	outputData.back()->SetData(device->GetLinacFlows()[0]->GetDynamicsAccel().GetData(),
//
device->GetOpticElementsSequence()->getParameter(device->GetOpticElementsSequence()->length()-1,
// "at")*1e-9, 0, 1,  0);
outputData.back()->addAccelElemetsDescription(device);
outputData.back()->SetAdditionalData(TimeOut, OutParams);

savetrack(outpuFilePaths[SynchrotronSolversBeam][0], device->GetOpticElementsSequence(), ind,
          outputData.back()->madXData);
progress = 1;*/
};

void SynchrotronSolver::DynamicsSimulationCenterMass(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{

    progress = 0;

    auto pt = readJSONFile(inputFileNames[SynchrotronSolversNameCenterOfMass][0]);

    if (!pt)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto errors = read_errors(*pt);

    if (!errors)
    {
        errorMsg = "Unable to read tolerances file";
        return;
    }

    auto seq_err = std::make_shared<OpticElementsSequence>(*device->GetOpticElementsSequence());

    seq_err->InsertErrors(*errors);

    seq_err->mtrans();

    auto result = simulationCenterMass<float>(device->GetLinacFlows()[0], seq_err,
                                              device->GetParameter("Number of cicles"), "ALL", 0);

    outputData.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(SynchrotronDynCMFlags, SynchrotronSolversNameCenterOfMass)));

    std::vector<std::vector<float>> OutParams(3);

    OutParams[0].resize(result[0].size());
    OutParams[1].resize(result[0].size());

    std::fill(OutParams[0].begin(), OutParams[0].end(), device->GetParameter("X aperture, m"));
    std::fill(OutParams[1].begin(), OutParams[1].end(), device->GetParameter("Y aperture, m"));

    outputData.back()->addAccelElemetsDescription(device);

    outputData.back()->SetAdditionalData(result[0], OutParams);

    outputData.back()->YData.resize(1);
    outputData.back()->YData[0] = result;

    progress = 1;
};

void SynchrotronSolver::SynchrotronSolversMagnetShuffle(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    std::vector<std::string>& status)
{
    outputData.push_back(std::shared_ptr<SimulationDataAccel>(
        new SimulationDataAccel(SynchrotronShuffleFlags, SynchrotronSolversMagnetShuffleName)));

    auto arr = device->GetOpticElementsSequence()->get_types_array("SBEND");

    auto n_iter =
        int(solverParameters[SynchrotronSolversMagnetShuffleName]->find("Number of iterations"));

    std::map<std::string, double> actual_vals;

    std::ifstream input_vals(inputFileNames[SynchrotronSolversMagnetShuffleName][0]);

    while (input_vals)
    {
        std::string tmp;
        double      tmp1;

        input_vals >> tmp;
        input_vals >> tmp1;
        actual_vals[tmp] = tmp1;
    }

    outputData.back()->XData.resize(2);

    for (auto& v : outputData.back()->XData)
    {
        v.resize(n_iter);
    }
    double r_min = 1e308;

    auto best_shuffle = arr;

    for (int i = 0; i < n_iter; i++)
    {
        std::random_shuffle(arr.begin(), arr.end());

        auto seq_new = std::make_shared<OpticElementsSequence>(*device->GetOpticElementsSequence());

        seq_new->insert_shuffle("SBEND", arr, actual_vals);
        seq_new->mtrans();

        auto result_orbit = calc_closed_orbit(
            seq_new, inputFileNames[SynchrotronSolversNameClosedOrbit][1],
            device->GetParameter("X aperture, m"), device->GetParameter("Y aperture, m"));

        device->GetLinacFlows()[0]->SetMassCenterVector(result_orbit.first.first, false);
        device->GetLinacFlows()[0]->SetMassCenterVector(result_orbit.first.second, true);

        auto result =
            simulationCenterMass<float>(device->GetLinacFlows()[0], seq_new,
                                        device->GetParameter("Number of cicles"), "ALL", 1);

        outputData.back()->XData[0][i] =
            std::abs(*std::max_element(result[1].begin(), result[1].end(), [](auto v1, auto v2) {
                return std::abs(v1) < std::abs(v2);
            }));
        outputData.back()->XData[1][i] =
            std::abs(*std::max_element(result[2].begin(), result[2].end(), [](auto v1, auto v2) {
                return std::abs(v1) < std::abs(v2);
            }));

        auto r = std::max(outputData.back()->XData[0][i], outputData.back()->XData[1][i]);

        if (r < r_min)
        {
            r_min        = r;
            best_shuffle = arr;
        }
        progress = double(i) / n_iter;
    }
    progress = 1.0;
    std::ofstream output_vals(outpuFilePaths[SynchrotronSolversMagnetShuffleName][0]);

    for (const auto& v : best_shuffle)
    {
        output_vals << v << "\t";
    }
}

void SynchrotronSolver::OrbitCalculation(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg){
    /* outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
         SynchrotronOrbitFlags, SynchrotronSolversOrbitCalc, 4, {1, 1, 1, 1})));

     std::vector<arma::mat> xtwiss;
     std::vector<arma::mat> ytwiss;

     std::vector<double> mu_y;
     std::vector<double> mu_x;

     std::vector<double> S;

     DynamicsSimulationTwissF(device, xtwiss, ytwiss, mu_x, mu_y, S);
     std::vector<arma::mat> x;
     std::vector<arma::mat> y;
     InitOrbit(device, x, y);

     std::vector<int> nump = device->GetOpticElementsSequence()->findType("MONITOR");

     std::vector<double> xn(nump.size());
     std::vector<double> yn(nump.size());

     for (int i = 0; i < nump.size(); i++)
     {
         xn[i] = x[i](0, 0) / (sqrt(xtwiss[i](0, 0)));
         yn[i] = y[i](0, 0) / (sqrt(ytwiss[i](0, 0)));
     }
     int                 n = device->GetOpticElementsSequence()->length();
     std::vector<double> sxn(n);
     std::vector<double> syn(n);
     std::vector<double> s(n);

     double L  = device->GetOpticElementsSequence()->GetL();
     double Qx = device->GetParameter("Frequency of betatron oscillations");
     for (int i = 0; i < n; i++)
     {
         s[i]   = S[i];
         sxn[i] = s[i] + mu_x[i] * L / (2 * commtools::PI() * Qx);
         syn[i] = s[i] + mu_y[i] * L / (2 * commtools::PI() * Qx);
     }
     std::vector<double> ax(nump.size() - 1);
     std::vector<double> phix(nump.size() - 1);
     std::vector<double> ay(nump.size() - 1);
     std::vector<double> phiy(nump.size() - 1);

     for (int i = 0; i < nump.size() - 1; i++)
     {
         eqsolve(ax[i], phix[i], xn[i], xn[i + 1], 2 * commtools::PI() * Qx * sxn[nump[i]] / L,
                 2 * commtools::PI() * Qx * sxn[nump[i + 1]] / L);
         eqsolve(ay[i], phiy[i], yn[i], yn[i + 1], 2 * commtools::PI() * Qx * syn[nump[i]] / L,
                 2 * commtools::PI() * Qx * syn[nump[i + 1]] / L);
     }
     std::vector<double> sxnp(nump.size());
     std::vector<double> synp(nump.size());
     for (int i = 0; i < nump.size(); i++)
     {
         sxnp[i] = sxn[nump[i]];
         synp[i] = syn[nump[i]];
     }
     std::vector<double> xx(n);
     std::vector<double> yy(n);

     std::vector<int> ix;
     std::vector<int> iy;

     for (int i = nump[0]; i < nump.back(); i++)
     {
         ix.clear();
         iy.clear();

         // int i = nump[k];
         tabl_ind3(ix, sxnp, sxn[i]);
         if (ix.size())
         {
             double x12 = ax[ix[0]] * std::sin(2 * commtools::PI() * Qx * sxn[i] / L + phix[ix[0]]);
             double x23 = ax[ix[1]] * std::sin(2 * commtools::PI() * Qx * sxn[i] / L + phix[ix[1]]);
             xx[i]      = calcx(x12, x23, sxn[i], sxnp[ix[0]], sxnp[ix[1]], sxnp[ix[2]]) *
                     sqrt(xtwiss[i](0, 0));
         }
         // outputData.back()->addData(0, S[i], xx[i]);

         tabl_ind3(iy, synp, syn[i]);
         if (iy.size())
         {
             double y12 = ay[iy[0]] * std::sin(2 * commtools::PI() * Qx * syn[i] / L + phiy[iy[0]]);
             double y23 = ay[iy[1]] * std::sin(2 * commtools::PI() * Qx * syn[i] / L + phiy[iy[1]]);
             yy[i]      = calcx(y12, y23, syn[i], synp[ix[0]], synp[iy[1]], synp[iy[2]]) *
                     sqrt(ytwiss[i](0, 0));
         }
         //	outputData.back()->addData(1, S[i], yy[i]);
     }

     std::vector<double> initialPos;
     device->GetLinacFlows()[0]->GetMassCenterVector(initialPos);

     arma::mat x1 = {initialPos[0], initialPos[1], 0};
     x1           = arma::trans(x1);
     arma::mat y1 = {initialPos[2], initialPos[3], 0};
     y1           = arma::trans(y1);
     int   nrb    = 10;
     int   nqd    = 10;
     float ss     = 0;

     std::vector<std::vector<float>> OutParams(3);
     std::vector<float>              TimeOut;

     outputData.back()->addData(0, ss, float(x1(0, 0)));
     outputData.back()->addData(1, ss, float(y1(0, 0)));

     for (int i = 0; i < device->GetOpticElementsSequence()->length(); i++)
     {
         if (device->GetOpticElementsSequence()->GetType(i) == "RBEND")
         {
             arma::mat Msh(3, 3);
             arma::mat Msv(3, 3);
             arma::mat Me1h(3, 3);
             arma::mat Me2h(3, 3);
             arma::mat Me1v(3, 3);
             arma::mat Me2v(3, 3);

             double L = device->GetOpticElementsSequence()->getParameter(i, "L") / nrb;
             double R = device->GetOpticElementsSequence()->getParameter(i, "L") /
                        device->GetOpticElementsSequence()->getParameter(i, "ANGLE");
             double eta  = L / R;
             double eps1 = device->GetOpticElementsSequence()->getParameter(i, "E1");
             double eps2 = device->GetOpticElementsSequence()->getParameter(i, "E2");

             double pp = L / eta;

             Msh(0, 0) = std::cos(eta);
             Msh(0, 1) = pp * std::sin(eta);
             Msh(0, 2) = pp * (1 - std::cos(eta));
             Msh(1, 0) = -std::sin(eta) / pp;
             Msh(1, 1) = std::cos(eta);
             Msh(1, 2) = std::sin(eta);
             Msh(2, 0) = 0;
             Msh(2, 1) = 0;
             Msh(2, 2) = 1;

             Msv(0, 0) = 1;
             Msv(0, 1) = pp * eta;
             Msv(0, 2) = 0;
             Msv(1, 0) = 0;
             Msv(1, 1) = 1;
             Msv(1, 2) = 0;
             Msv(2, 0) = 0;
             Msv(2, 1) = 0;
             Msv(2, 2) = 1;

             eta = device->GetOpticElementsSequence()->getParameter(i, "ANGLE");
             pp  = device->GetOpticElementsSequence()->getParameter(i, "L") / eta;

             Me1h(0, 0) = 1;
             Me1h(0, 1) = 0;
             Me1h(0, 2) = 0;
             Me1h(1, 0) = tan(eta / 2 + eps1) / pp;
             Me1h(1, 1) = 1;
             Me1h(1, 2) = 0;
             Me1h(2, 0) = 0;
             Me1h(2, 1) = 0;
             Me1h(2, 2) = 1;

             Me2h(0, 0) = 1;
             Me2h(0, 1) = 0;
             Me2h(0, 2) = 0;
             Me2h(1, 0) = tan(eta / 2 + eps2) / pp;
             Me2h(1, 1) = 1;
             Me2h(1, 2) = 0;
             Me2h(2, 0) = 0;
             Me2h(2, 1) = 0;
             Me2h(2, 2) = 1;

             Me1v(0, 0) = 1;
             Me1v(0, 1) = 0;
             Me1v(0, 2) = 0;
             Me1v(1, 0) = -tan(eta / 2 + eps1) / pp;
             Me1v(1, 1) = 1;
             Me1v(1, 2) = 0;
             Me1v(2, 0) = 0;
             Me1v(2, 1) = 0;
             Me1v(2, 2) = 1;

             Me2v(0, 0) = 1;
             Me2v(0, 1) = 0;
             Me2v(0, 2) = 0;
             Me2v(1, 0) = -tan(eta / 2 + eps2) / pp;
             Me2v(1, 1) = 1;
             Me2v(1, 2) = 0;
             Me2v(2, 0) = 0;
             Me2v(2, 1) = 0;
             Me2v(2, 2) = 1;

             x1 = Msh * Me1h * x1;
             y1 = Msv * Me1v * y1;
             ss = ss + L;

             outputData.back()->addData(0, ss, float(x1(0, 0)));
             outputData.back()->addData(1, ss, float(y1(0, 0)));

             for (int j = 1; j < nrb - 1; j++)
             {
                 x1 = Msh * x1;
                 y1 = Msv * y1;
                 ss = ss + L;

                 outputData.back()->addData(0, ss, float(x1(0, 0)));
                 outputData.back()->addData(1, ss, float(y1(0, 0)));
             }

             x1 = Me2h * Msh * x1;
             y1 = Me2v * Msv * y1;
             ss = ss + L;

             outputData.back()->addData(0, ss, float(x1(0, 0)));
             outputData.back()->addData(1, ss, float(y1(0, 0)));
             OutParams[0].push_back(device->GetParameter("X aperture, m"));
             OutParams[1].push_back(device->GetParameter("Y aperture, m"));
             TimeOut.push_back(ss);
             continue;
         };
         if (device->GetOpticElementsSequence()->GetType(i) == "QUADRUPOLE")
         {
             double    Kx   = std::abs(device->GetOpticElementsSequence()->getParameter(i, "k1"));
             double    Ky   = std::abs(device->GetOpticElementsSequence()->getParameter(i, "k1"));
             double    Lm   = device->GetOpticElementsSequence()->getParameter(i, "L") / nqd;
             double    ksix = Lm * sqrt(Kx);
             double    ksiy = Lm * sqrt(Ky);
             arma::mat Tx   = arma::mat(3, 3, arma::fill::eye);
             arma::mat Ty   = arma::mat(3, 3, arma::fill::eye);
             if (device->GetOpticElementsSequence()->getParameter(i, "k1") > 0) // ������������
     �����
             {

                 Tx(0, 0) = std::cos(ksix);
                 Tx(0, 1) = 1.0 / sqrt(Kx) * std::sin(ksix);
                 Tx(0, 2) = 0;
                 Tx(1, 0) = -sqrt(Kx) * std::sin(ksix);
                 Tx(1, 1) = std::cos(ksix);
                 Tx(1, 2) = 0;
                 Tx(2, 0) = 0;
                 Tx(2, 1) = 0;
                 Tx(2, 2) = 1;

                 Ty(0, 0) = cosh(ksiy);
                 Ty(0, 1) = 1.0 / sqrt(Ky) * sinh(ksiy);
                 Ty(0, 2) = 0;
                 Ty(1, 0) = sqrt(Ky) * sinh(ksiy);
                 Ty(1, 1) = cosh(ksiy);
                 Ty(1, 2) = 0;
                 Ty(2, 0) = 0;
                 Ty(2, 1) = 0;
                 Ty(2, 2) = 1;
             }
             else // �������������� �����
             {
                 Tx(0, 0) = cosh(ksix);
                 Tx(0, 1) = 1.0 / sqrt(Kx) * sinh(ksix);
                 Tx(0, 2) = 0;
                 Tx(1, 0) = sqrt(Kx) * sinh(ksix);
                 Tx(1, 1) = cosh(ksix);
                 Tx(1, 2) = 0;
                 Tx(2, 0) = 0;
                 Tx(2, 1) = 0;
                 Tx(2, 2) = 1;

                 Ty(0, 0) = std::cos(ksiy);
                 Ty(0, 1) = 1.0 / sqrt(Ky) * std::sin(ksiy);
                 Ty(0, 2) = 0;
                 Ty(1, 0) = -sqrt(Ky) * std::sin(ksiy);
                 Ty(1, 1) = std::cos(ksiy);
                 Ty(1, 2) = 0;
                 Ty(2, 0) = 0;
                 Ty(2, 1) = 0;
                 Ty(2, 2) = 1;
             };
             for (int j = 0; j < nqd; j++)
             {
                 x1 = Tx * x1;
                 y1 = Ty * y1;
                 ss = ss + Lm;
                 outputData.back()->addData(0, ss, float(x1(0, 0)));
                 outputData.back()->addData(1, ss, float(y1(0, 0)));
             }
             OutParams[0].push_back(device->GetParameter("X aperture, m"));
             OutParams[1].push_back(device->GetParameter("Y aperture, m"));
             TimeOut.push_back(ss);
             continue;
         }
         x1 = device->GetOpticElementsSequence()->Tx[i] * x1;
         y1 = device->GetOpticElementsSequence()->Tx[i] * y1;
         ss = ss + device->GetOpticElementsSequence()->getParameter(i, "L");
         outputData.back()->addData(0, ss, float(x1(0, 0)));
         outputData.back()->addData(1, ss, float(y1(0, 0)));

         OutParams[0].push_back(device->GetParameter("X aperture, m"));
         OutParams[1].push_back(device->GetParameter("Y aperture, m"));
         TimeOut.push_back(ss);
     }
     outputData.back()->addAccelElemetsDescription(device);
     outputData.back()->SetAdditionalData(TimeOut, OutParams);*/
};

void SynchrotronSolver::AcceptanceCalculation(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{
    outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
        SynchrotronAccFlags, SynchrotronSolversAcceptance, 4, {0, 0, 0, 0})));
    // outputData.push_back(std::shared_ptr<SimulationDataAccel>(new
    // SimulationDataAccel(SynchrotronDynFlags, solversNames[2], 6, { 1,1,1,1,1,1 })));

    int ncircle = device->GetParameter("Number of cicles");
    device->GetLinacFlows()[0]->GenerateParticlesAccelTest(
        solverParameters[SynchrotronSolversAcceptance]->find("Number of particles"),
        device->GetParameter("X aperture, m"), device->GetParameter("Y aperture, m"),
        solverParameters[SynchrotronSolversAcceptance]->find("dX/dZ search range"),
        solverParameters[SynchrotronSolversAcceptance]->find("dY/dZ search range"));

    double rx     = device->GetParameter("X aperture, m");
    double ry     = device->GetParameter("Y aperture, m");
    double rel_ch = solverParameters[SynchrotronSolversAcceptance]->find("Relative channel radius");

    device->GetLinacFlows()[0]->GetDynamicsAccelTest().Init();

    int ns = 5000;

    double atC = 0;
    size_t k   = 0;

    auto m_indexes = getIndexes(device->GetOpticElementsSequence(), "ALL");

    std::vector<float> xStart;
    std::vector<float> dxStart;
    std::vector<float> yStart;
    std::vector<float> dyStart;

    for (auto& x0 : device->GetLinacFlows()[0]->GetDynamicsAccelTest().positionStart)
    {
        progress =
            float(k) / (device->GetLinacFlows()[0]->GetDynamicsAccelTest().positionStart.size());
        std::vector<std::vector<float>> result =
            simulationTrace<float>(x0, device->GetOpticElementsSequence(), ncircle, m_indexes, 1);

        bool flagAdd = true;

        for (int i = 0; i < result[0].size(); i++)
        {
            if ((result[1][i] * result[1][i]) / (rx * rx) +
                    (result[2][i] * result[2][i]) / (ry * ry) >
                rel_ch)
            {
                flagAdd = false;
                break;
            }
            else
            {
                int tt = 0;
            }
        }
        if (flagAdd)
        {
            xStart.push_back(x0(0));
            dxStart.push_back(x0(1));
            yStart.push_back(x0(2));
            dyStart.push_back(x0(3));
        }
    }

    outputData.back()->addData(0, 0, xStart);
    outputData.back()->addData(1, 0, dxStart);
    outputData.back()->addData(2, 0, yStart);
    outputData.back()->addData(3, 0, dyStart);

    /* for (int j = 0; j < ncircle; j++)
     {
         double at = 0;
         for (int i = 0; i < device->GetOpticElementsSequence()->length(); i++)
         {
             updateBeamPositions(device->GetOpticElementsSequence()->Tx[i],
                                 device->GetOpticElementsSequence()->Ty[i],
                                 device->GetLinacFlows()[0]->GetDynamicsAccelTest());
             at = atC + device->GetOpticElementsSequence()->getParameter(i, "at") +
                  device->GetOpticElementsSequence()->getParameter(i, "L");

             progress = double(j * device->GetOpticElementsSequence()->length() + i) /
                        double(ncircle * device->GetOpticElementsSequence()->length());
             device->GetLinacFlows()[0]->GetDynamicsAccelTest().CalculateTransmission(
                 device->GetParameter("X aperture, m"), device->GetParameter("Y aperture, m"),
                 solverParameters[SynchrotronSolversAcceptance]->find("Relative channel radius"));
             //
     outputData.back()->SetData(device->GetLinacFlows()[0]->GetDynamicsAccelTest().GetData(),
              at*1e-9, 0, 1,
              0);
         };
         atC = at;
     }
     device->GetLinacFlows()[0]->GetDynamicsAccelTest().removeParticle();

*/
};