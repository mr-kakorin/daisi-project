#include <armadillo>

#include "Results.h"
#include "SynchrotronDevice.h"
#include "SynchrotronSolver.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "SynchrotronToolsSim.h"
#include "Tools.h"

#include "../base/AccelFlow.h"

#include <common_tools/constants.h>

void invSVD(const arma::mat& R, const arma::mat& s, arma::mat& S)
{
    int p = R.n_rows;
    int k = R.n_cols;
    S     = (arma::mat(p, k, arma::fill::zeros));

    int l = std::min(p, k);
    for (int i = 0; i < l; i++)
    {
        if (s(i) != 0)
            S(i, i) = 1.0 / s(i);
    };
};
void SynchrotronSolver::Mikado(int flag, std::shared_ptr<SynchrotronDevice>& device, arma::mat& R,
                               std::vector<int>& k1, int N, arma::mat& dx0, const arma::mat& x0)
{
    std::vector<int> numk = device->GetOpticElementsSequence()->findType("KICKER");
    std::vector<int> nump = device->GetOpticElementsSequence()->findType("MONITOR");
    int              k    = numk.size(); // ����� �����������
    int              p    = nump.size(); // ����� �������

    arma::mat                              R1;
    std::shared_ptr<OpticElementsSequence> seqj =
        std::shared_ptr<OpticElementsSequence>(new OpticElementsSequence());
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::mat S;
    arma::vec vdx(numk.size());
    for (int k = 0; k < N; k++)
    {
        arma::mat Rt;
        for (int i = 0; i < numk.size(); i++)
        {
            Rt = arma::join_rows(R1, R.col(i));
            arma::svd(U, s, V, Rt);
            invSVD(Rt, s, S);
            arma::vec dI = V * arma::trans(S) * trans(U) * trans(dx0.row(0));
            *seqj        = *device->GetOpticElementsSequence();

            std::vector<int> kt = k1;
            kt.push_back(numk[i]);
            int    j;
            double hk, vk;
            for (j = 0; j < dI.n_rows; j++)
            {
                switch (flag)
                {
                case 0:
                    hk = seqj->getParameter(kt[j], "HKICK");
                    seqj->setParameter(kt[j], "HKICK", hk - dI(j));
                    break;
                case 1:
                    vk = seqj->getParameter(kt[j], "VKICK");
                    seqj->setParameter(kt[j], "VKICK", vk - dI(j));
                    break;
                }
            }
            seqj->mtrans();

            arma::mat dx = (arma::mat(3, p, arma::fill::zeros));

            for (int ii    = 0; ii < p; ii++)
                dx.col(ii) = seqj->sTx[nump[ii]] * x0;

            vdx(i) = arma::sum(dx.row(0) % dx.row(0));
        };

        int jx = vdx.index_min();
        R1     = arma::join_rows(R1, R.col(jx));
        k1.push_back(numk[jx]);
        R.shed_col(jx);
        numk.erase(numk.begin() + jx);
    };
    R = R1;
};

void SynchrotronSolver::CorrMatrixCalc(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{
    std::vector<double> initialPos;
    device->GetLinacFlows()[0]->GetMassCenterVector(initialPos);

    arma::mat x0 = {{initialPos[0], initialPos[1], 0}};
    x0           = arma::trans(x0);
    arma::mat y0 = {{initialPos[2], initialPos[3], 0}};
    y0           = arma::trans(y0);

    arma::mat Rx;
    arma::mat Ry;

    CorrMatrix(arma::vec(initialPos), device->GetOpticElementsSequence(), Rx, Ry,
               std::vector<int>());
    outpuFilePaths[SynchrotronSolversCorrMatrix][0];

    arma::mat R = arma::join_rows(Rx, Ry);
    R.save(outpuFilePaths[SynchrotronSolversCorrMatrix][0], arma::arma_ascii);
};
void SynchrotronSolver::OrbitCorrectionSVD(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{
    outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
        SynchrotronCorrFlags, SynchrotronSolversSVD, 6, {1, 1, 1, 1, 1, 1})));

    std::vector<int> numk = device->GetOpticElementsSequence()->findType("KICKER");
    std::vector<int> nump = device->GetOpticElementsSequence()->findType("MONITOR");
    int              k    = numk.size(); // ����� �����������
    int              p    = nump.size(); // ����� �������

    if (!k || !p)
    {
        errorMsg = "There are no pick-up electrodes or correctors in accelerator optic";
        return;
    };

    std::vector<double> initialPos;
    device->GetLinacFlows()[0]->GetMassCenterVector(initialPos);

    arma::mat x0 = {{initialPos[0], initialPos[1], 0}};
    x0           = arma::trans(x0);
    arma::mat y0 = {{initialPos[2], initialPos[3], 0}};
    y0           = arma::trans(y0);

    arma::mat dx0 = (arma::mat(3, p, arma::fill::zeros));
    arma::mat dy0 = (arma::mat(3, p, arma::fill::zeros));

    std::vector<std::vector<double>> res = simulationTrace<double>(
        arma::vec(initialPos), device->GetOpticElementsSequence(), 1, nump, 1);

    dx0.row(0) = arma::Row<double>(res[1]);
    dy0.row(0) = arma::Row<double>(res[2]);

    arma::mat Rx;
    arma::mat Ry;
    // CorrMatrix(device, Rx, Ry, x0, y0);

    arma::mat Rt;
    Rt.load(inputFileNames[SynchrotronSolversSVD][0]);

    int ncol = Rt.n_cols / 2;
    Rx       = Rt.cols(0, ncol - 1);
    Ry       = Rt.cols(ncol, ncol * 2 - 1);

    arma::mat Ux;
    arma::vec sx;
    arma::mat Vx;

    arma::mat Uy;
    arma::vec sy;
    arma::mat Vy;

    arma::mat Sx1;
    arma::mat Sy1;

    arma::svd(Ux, sx, Vx, Rx);
    arma::svd(Uy, sy, Vy, Ry);

    invSVD(Rx, sx, Sx1);
    invSVD(Ry, sy, Sy1);

    FILE* fp;
    fp = fopen(outpuFilePaths[SynchrotronSolversSVD][0].c_str(), "w");
    fprintf(fp, "Singular values for X correction matrix:\n");

    for (int i = 0; i < sx.n_rows; i++)
        fprintf(fp, "%lf\t", sx(i));

    fprintf(fp, "\nSingular values for Y correction matrix:\n");

    for (int i = 0; i < sx.n_rows; i++)
        fprintf(fp, "%lf\t", sy(i));

    fclose(fp);

    std::shared_ptr<OpticElementsSequence> seqj =
        std::shared_ptr<OpticElementsSequence>(new OpticElementsSequence());
    arma::vec dIx = Vx * arma::trans(Sx1) * arma::trans(Ux) * arma::trans(dx0.row(0));
    arma::vec dIy = Vy * arma::trans(Sy1) * arma::trans(Uy) * arma::trans(dy0.row(0));
    dIx.save("dIx.txt", arma::arma_ascii);
    dx0.save("dx0.txt", arma::arma_ascii);
    Vx.save("Vx.txt", arma::arma_ascii);
    Sx1.save("Sx1.txt", arma::arma_ascii);
    Ux.save("Ux.txt", arma::arma_ascii);
    Rx.save("Rx.txt", arma::arma_ascii);

    //   *seqj = *device->GetOpticElementsSequence();
    for (int j = 0; j < numk.size(); j++)
    {
        double hk = device->GetOpticElementsSequence()->getParameter(numk[j], "HKICK");
        double vk = device->GetOpticElementsSequence()->getParameter(numk[j], "VKICK");
        device->GetOpticElementsSequence()->setParameter(numk[j], "HKICK", hk - 1e-4 * dIx(j));
        device->GetOpticElementsSequence()->setParameter(numk[j], "VKICK", vk - 1e-4 * dIy(j));
    }

    device->GetOpticElementsSequence()->mtrans();

    arma::mat dx = (arma::mat(3, p, arma::fill::zeros));
    arma::mat dy = (arma::mat(3, p, arma::fill::zeros));

    std::vector<std::vector<double>> res_new = simulationTrace<double>(
        arma::vec(initialPos), device->GetOpticElementsSequence(), 1, nump, 1);

    dx.row(0) = arma::Row<double>(res_new[1]);
    dy.row(0) = arma::Row<double>(res_new[2]);

    outputData.back()->SetDataNucl(dx, dy, dx0, dy0, dIx, dIy, sx, sy);
};
void SynchrotronSolver::OrbitCorrectionMikado(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{
    outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
        SynchrotronCorrFlags, SynchrotronSolversMikado, 6, {1, 1, 1, 1, 1, 1})));

    std::vector<int>    numk = device->GetOpticElementsSequence()->findType("KICKER");
    std::vector<int>    nump = device->GetOpticElementsSequence()->findType("MONITOR");
    int                 k    = numk.size(); // ����� �����������
    int                 p    = nump.size(); // ����� �������
    std::vector<double> initialPos;
    device->GetLinacFlows()[0]->GetMassCenterVector(initialPos);
    // arma::mat x0 = { {3e-3, -2e-5, 0} };
    arma::mat x0 = {{initialPos[0], initialPos[1], 0}};

    x0 = arma::trans(x0);
    // arma::mat y0 = { { -3e-3, -2e-5, 0} };
    arma::mat y0 = {{initialPos[2], initialPos[3], 0}};

    y0 = arma::trans(y0);

    arma::mat dx0 = (arma::mat(3, p, arma::fill::zeros));
    arma::mat dy0 = (arma::mat(3, p, arma::fill::zeros));
    for (int i = 0; i < p; i++)
    {
        dx0.col(i) = device->GetOpticElementsSequence()->sTx[nump[i]] * x0;
        dy0.col(i) = device->GetOpticElementsSequence()->sTy[nump[i]] * y0;
    }
    arma::mat Rx;
    arma::mat Ry;
    // CorrMatrix(device, Rx, Ry, x0, y0);

    arma::mat Rt;
    Rt.load(inputFileNames[SynchrotronSolversMikado][0]);

    int ncol = Rt.n_cols / 2;
    Rx       = Rt.cols(0, ncol - 1);
    Ry       = Rt.cols(ncol, ncol * 2 - 1);

    int Nkx = int(solverParameters[SynchrotronSolversMikado]->find("Number of correctors (X)"));
    int Nky = int(solverParameters[SynchrotronSolversMikado]->find("Number of correctors (Y)"));

    std::vector<int> numkx;
    std::vector<int> numky;

    Mikado(0, device, Rx, numkx, Nkx, dx0, x0);
    Mikado(1, device, Ry, numky, Nky, dy0, y0);

    arma::mat Ux;
    arma::vec sx;
    arma::mat Vx;

    arma::mat Uy;
    arma::vec sy;
    arma::mat Vy;

    arma::mat Sx1;
    arma::mat Sy1;

    arma::svd(Ux, sx, Vx, Rx);
    arma::svd(Uy, sy, Vy, Ry);

    invSVD(Rx, sx, Sx1);
    invSVD(Ry, sy, Sy1);

    FILE* fp = fopen(outpuFilePaths[SynchrotronSolversMikado][0].c_str(), "w");
    fprintf(fp, "Singular values for X correction matrix:\n");

    for (int i = 0; i < sx.n_rows; i++)
        fprintf(fp, "%lf\t", sx(i));

    fprintf(fp, "\nSingular values for Y correction matrix:\n");

    for (int i = 0; i < sx.n_rows; i++)
        fprintf(fp, "%lf\t", sy(i));

    fclose(fp);

    std::shared_ptr<OpticElementsSequence> seqj =
        std::shared_ptr<OpticElementsSequence>(new OpticElementsSequence());
    arma::vec dIx = Vx * arma::trans(Sx1) * arma::trans(Ux) * arma::trans(dx0.row(0));
    arma::vec dIy = Vy * arma::trans(Sy1) * arma::trans(Uy) * arma::trans(dy0.row(0));
    // dIx.save("dIx.txt",arma::arma_ascii);
    *seqj = *device->GetOpticElementsSequence();
    for (int j = 0; j < dIx.n_rows; j++)
    {
        double hk = seqj->getParameter(numkx[j], "HKICK");
        seqj->setParameter(numkx[j], "HKICK", hk - dIx(j));
    }
    for (int j = 0; j < dIy.n_rows; j++)
    {
        double vk = seqj->getParameter(numky[j], "VKICK");
        seqj->setParameter(numky[j], "VKICK", vk - dIy(j));
    }

    seqj->mtrans();

    arma::mat dx = (arma::mat(3, p, arma::fill::zeros));
    arma::mat dy = (arma::mat(3, p, arma::fill::zeros));

    for (int i = 0; i < p; i++)
    {
        dx.col(i) = seqj->sTx[nump[i]] * x0;
        dy.col(i) = seqj->sTy[nump[i]] * y0;
    }
    outputData.back()->SetDataNucl(dx, dy, dx0, dy0, dIx, dIy, sx, sy);
};
