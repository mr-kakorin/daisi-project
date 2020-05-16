#include "SynchrotronToolsSim.h"
#include "../base/AccelFlow.h"
#include "SynchrotronDevice.h"
#include "SynchrotronTools.h"

#include <common_tools/constants.h>
#include <notk/controller.hpp>

#include <algorithm>
#include <memory>
#include <vector>

int matrType = 0;

std::vector<int> getIndexes(const std::shared_ptr<OpticElementsSequence>& sequence,
                            const std::string&                            regime)
{
    std::vector<int> m_indexes;
    if (regime == "ALL")
    {
        m_indexes.resize(sequence->length());
        for (size_t i = 0; i < m_indexes.size(); i++)
        {
            m_indexes[i] = i;
        }
    }
    if (regime == "MONITOR")
    {
        m_indexes = sequence->findType("MONITOR");
    }
    if (regime == "FIRST&LAST")
    {
        m_indexes.resize(2);
        m_indexes[0] = 1;
        m_indexes[1] = sequence->length() - 1;
    }
    return m_indexes;
}

template std::vector<std::vector<float>> DynamicsSimulationTwissF<float>(
    const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle);

template std::vector<std::vector<double>> DynamicsSimulationTwissF<double>(
    const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle);

template <class T>
std::vector<std::vector<T>>
DynamicsSimulationTwissF(const arma::vec&                              x0,
                         const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle)
{

    T lasty = 0;
    T lastx = 0;

    arma::mat                   Mx(2, 2, arma::fill::eye);
    arma::mat                   My(2, 2, arma::fill::eye);
    std::vector<std::vector<T>> dyn(9);
    for (auto& d : dyn)
    {
        d.resize(ncircle * sequence->length());
    }
    double atC = 0;

    if (matrType == 1 || matrType == 0)
    {
        arma::vec xs = {x0[0], x0[1], x0[2]};
        arma::vec ys = {x0[3], x0[4], x0[5]};

        double bx1 = xs(0, 0);
        double by1 = ys(0, 0);

        size_t k = 0;

        for (int j = 0; j < ncircle; j++)
        {
            for (int i = 0; i < sequence->length(); i++)
            {
                Mx = sequence->Tx[i].submat(0, 0, 1, 1);
                My = sequence->Ty[i].submat(0, 0, 1, 1);

                xs = twiss(Mx, xs);
                ys = twiss(My, ys);

                double bx2  = xs(0);
                double by2  = ys(0);
                double mx12 = sequence->Tx[i](0, 1);
                double my12 = sequence->Ty[i](0, 1);
                double mu_x = (lastx + std::asin(mx12 / sqrt(bx1 * bx2)) / (2 * commtools::PI()));
                double mu_y = (lasty + std::asin(my12 / sqrt(by1 * by2)) / (2 * commtools::PI()));
                lasty       = mu_x;
                lastx       = mu_y;
                bx1         = bx2;
                by1         = by2;

                atC = sequence->getParameter(i, "at") + sequence->getParameter(i, "L") +
                      j * sequence->GetL();

                dyn[0][k] = atC;

                dyn[1][k] = xs(0);
                dyn[2][k] = xs(1);
                dyn[3][k] = xs(2);

                dyn[4][k] = ys(0);
                dyn[5][k] = ys(1);
                dyn[6][k] = ys(2);

                dyn[7][k] = mu_x;
                dyn[8][k] = mu_y;

                k++;
            };
        }
    }
    return dyn;
}

void InitOrbit(std::shared_ptr<SynchrotronDevice>& device, std::vector<arma::mat>& x,
               std::vector<arma::mat>& y)
{
    std::vector<double> initialPos;
    device->GetLinacFlows()[0]->GetMassCenterVector(initialPos);

    arma::mat x1 = {initialPos[0], initialPos[1], 0};
    x1           = arma::trans(x1);
    arma::mat y1 = {initialPos[2], initialPos[3], 0};
    y1           = arma::trans(y1);

    std::vector<int> nump = device->GetOpticElementsSequence()->findType("MONITOR");
    x.resize(nump.size());
    y.resize(nump.size());

    for (int i = 0; i < nump.size(); i++)
    {
        x[i] = device->GetOpticElementsSequence()->sTx[nump[i]] * x1;
        y[i] = device->GetOpticElementsSequence()->sTy[nump[i]] * y1;
    }
}

void CorrMatrix(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& seq,
                arma::mat& Rx, arma::mat& Ry, std::vector<int> nump)
{
    std::vector<int> numk = seq->findType("KICKER");
    if (nump.size() == 0)
    {
        nump = getIndexes(seq, "MONITOR");
    }
    int k = numk.size(); // ����� �����������
    int p = nump.size(); // ����� �������

    double alpha = 1e-4;

    Rx = (arma::mat(p, k, arma::fill::zeros));
    Ry = (arma::mat(p, k, arma::fill::zeros));

    std::shared_ptr<OpticElementsSequence> seqj(std::make_shared<OpticElementsSequence>());

    for (int j = 0; j < numk.size(); j++)
    {
        *seqj     = *seq;
        double hk = seq->getParameter(numk[j], "HKICK");
        double vk = seq->getParameter(numk[j], "VKICK");
        seqj->setParameter(numk[j], "HKICK", hk + alpha);
        seqj->setParameter(numk[j], "VKICK", vk + alpha);
        seqj->mtrans();

        std::vector<std::vector<double>> res = simulationTrace<double>(x0, seqj, 1, nump, 1);

        Rx.col(j) = arma::vec(res[1]);
        Ry.col(j) = arma::vec(res[2]);
    }
}

template std::vector<std::vector<float>>
simulationTrace<float>(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence,
                       const int ncircle, std::vector<int> m_indexes, const int saveType);

template std::vector<std::vector<double>>
simulationTrace<double>(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence,
                        const int ncircle, std::vector<int> m_indexes, const int saveType);

template <class Tout>
std::vector<std::vector<Tout>>
simulationTrace(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence,
                const int ncircle, std::vector<int> m_indexes, const int saveType)
{
    std::vector<std::vector<Tout>> dyn;
    if (saveType == 0)
    {
        dyn.resize(6);
    }
    if (saveType == 1)
    {
        dyn.resize(3);
    }
    for (auto& d : dyn)
    {
        d.resize(ncircle * m_indexes.size());
    }
    double atC = 0;

    if (matrType == 1)
    {
        arma::vec x0L = {x0[0], x0[1], x0[4]};
        arma::vec y0L = {x0[2], x0[3], x0[4]};
        arma::vec xs;
        arma::vec ys;
        size_t    k = 0;
        for (int j = 0; j < ncircle; j++)
        {
            for (int i = 0; i < m_indexes.size(); i++)
            {
                xs = sequence->sTx[m_indexes[i]] * x0L + sequence->sBx[m_indexes[i]];
                ys = sequence->sTy[m_indexes[i]] * y0L + sequence->sBy[m_indexes[i]];

                atC = sequence->getParameter(i, "at") + sequence->getParameter(i, "L") +
                      j * sequence->GetL();

                dyn[0][k] = atC;

                if (saveType == 0)
                {
                    dyn[1][k] = xs(0);
                    dyn[2][k] = xs(1);
                    dyn[3][k] = ys(0);
                    dyn[4][k] = ys(1);
                }
                if (saveType == 1)
                {
                    dyn[1][k] = xs(0);
                    dyn[2][k] = ys(0);
                }
                k++;
            }
        }
    }
    if (matrType == 0)
    {
        arma::vec xs;
        size_t    k = 0;

        for (int j = 0; j < ncircle; j++)
        {
            for (int i = 0; i < m_indexes.size(); i++)
            {
                xs = sequence->sT[m_indexes[i]] * x0 + sequence->sB[m_indexes[i]];

                atC = sequence->getParameter(i, "at") + sequence->getParameter(i, "L") +
                      j * sequence->GetL();

                dyn[0][k] = atC;

                if (saveType == 0)
                {
                    dyn[1][k] = xs(0);
                    dyn[2][k] = xs(1);
                    dyn[3][k] = xs(2);
                    dyn[4][k] = xs(3);
                }
                if (saveType == 1)
                {
                    dyn[1][k] = xs(0);
                    dyn[2][k] = xs(2);
                }
                k++;
            }
        }
    }
    return dyn;
}

template std::vector<std::vector<float>>
simulationCenterMass<float>(const std::vector<double>&                    initialPos,
                            const std::shared_ptr<OpticElementsSequence>& sequence,
                            const int ncircle, const std::string& regime, const int saveType);

template std::vector<std::vector<double>>
simulationCenterMass<double>(const std::vector<double>&                    initialPos,
                             const std::shared_ptr<OpticElementsSequence>& sequence,
                             const int ncircle, const std::string& regime, const int saveType);

template <class Tout>
std::vector<std::vector<Tout>>
simulationCenterMass(const std::vector<double>&                    initialPos,
                     const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle,
                     const std::string& regime, const int saveType)
{
    std::vector<int> m_indexes = getIndexes(sequence, regime);
    return simulationTrace<Tout>(arma::vec(initialPos), sequence, ncircle, m_indexes, saveType);
}

template std::vector<std::vector<float>>
simulationCenterMass<float>(std::shared_ptr<SynchrotronFlow>&             flow,
                            const std::shared_ptr<OpticElementsSequence>& sequence,
                            const int ncircle, const std::string& regime, const int saveType);

template std::vector<std::vector<double>>
simulationCenterMass<double>(std::shared_ptr<SynchrotronFlow>&             flow,
                             const std::shared_ptr<OpticElementsSequence>& sequence,
                             const int ncircle, const std::string& regime, const int saveType);

template <class Tout>
std::vector<std::vector<Tout>>
simulationCenterMass(std::shared_ptr<SynchrotronFlow>&             flow,
                     const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle,
                     const std::string& regime, const int saveType)
{
    std::vector<double> initialPos;
    flow->GetMassCenterVector(initialPos);

    return simulationCenterMass<Tout>(initialPos, sequence, ncircle, regime, saveType);
}

template void simulationBeam<float>(std::vector<std::vector<float>>&              Xdata,
                                    std::vector<std::vector<std::vector<float>>>& Ydata,
                                    std::shared_ptr<SynchrotronFlow>&             flow,
                                    const std::shared_ptr<OpticElementsSequence>& sequence,
                                    const int ncircle, double& progress, bool& flagAbort,
                                    const std::string& regime, const int saveType,
                                    const int nParticles);

template void simulationBeam<double>(std::vector<std::vector<float>>&              Xdata,
                                     std::vector<std::vector<std::vector<float>>>& Ydata,
                                     std::shared_ptr<SynchrotronFlow>&             flow,
                                     const std::shared_ptr<OpticElementsSequence>& sequence,
                                     const int ncircle, double& progress, bool& flagAbort,
                                     const std::string& regime, const int saveType,
                                     const int nParticles);

template <class Tout>
void simulationBeam(std::vector<std::vector<float>>&              Xdata,
                    std::vector<std::vector<std::vector<float>>>& Ydata,
                    std::shared_ptr<SynchrotronFlow>&             flow,
                    const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle,
                    double& progress, bool& flagAbort, const std::string& regime,
                    const int saveType, const int nParticles)
{
    std::vector<int> m_indexes = getIndexes(sequence, regime);
    Ydata.resize(4);
    Xdata.resize(nParticles);
    std::for_each(Ydata.begin(), Ydata.end(), [&](auto& vec) { vec.resize(nParticles); });
    size_t k = 0;
    for (auto& x0 : flow->GetDynamicsAccel().positionStart)
    {
        progress = float(k) / (flow->GetDynamicsAccel().positionStart.size());
        std::vector<std::vector<float>> result =
            simulationTrace<float>(x0, sequence, ncircle, m_indexes, saveType);

        Xdata[k] = result[0];
        for (int s = 0; s < 4; s++)
            Ydata[s][k] = result[s + 1];
        k++;
        if (k == nParticles)
        {
            break;
        }
    }
}
std::vector<float> calculateTransmission(const std::vector<std::vector<std::vector<float>>>& Ydata,
                                         float rx, float ry, float ChannelRel)
{
    std::vector<float> result;
    std::vector<bool>  isInPerture(Ydata[0].size());
    std::fill(isInPerture.begin(), isInPerture.end(), true);

    for (size_t j = 0; j < Ydata[0][0].size(); j++)
    {
        double tr = 0;
        for (int i = 0; i < Ydata[0].size(); i++)
        {
            if (!isInPerture[i])
                continue;

            if ((Ydata[0][i][j] * Ydata[0][i][j]) / (rx * rx) +
                    (Ydata[2][i][j] * Ydata[2][i][j]) / (ry * ry) <
                ChannelRel)
                tr++;
            else
                isInPerture[i] = false;
        };
        result.push_back(tr / Ydata[0].size());
    }

    return result;
}

std::pair<std::pair<std::vector<double>, std::vector<double>>, std::string>
calc_closed_orbit(const std::shared_ptr<OpticElementsSequence>& sequence,
                  const std::string& notk_config, const double x_max, const double y_max)
{
    std::pair<std::pair<std::vector<double>, std::vector<double>>, std::string> ret_result;

	/* auto& errorMsg = ret_result.second;

    auto get_rel = [](const double v1, const double v2) -> double {
        if (v2 != 0)
        {
            return std::abs((v1 - v2) / v2);
        }
        else if (v1 != 0)
        {
            return std::abs((v1 - v2) / v1);
        }
        return std::abs(v1 - v2);
    };

    auto fitness = [&](const std::vector<double>& cm, const bool is_x) -> double {

        std::vector<double> cm_cur(5);
        if (is_x)
        {
            cm_cur[0] = cm[0];
            cm_cur[1] = cm[1];
        }
        else
        {
            cm_cur[2] = cm[0];
            cm_cur[3] = cm[1];
        }
        auto result = simulationCenterMass<double>(cm_cur, sequence, 1, "FIRST&LAST", 0);

        if (is_x)
        {
            return get_rel(result[1].back(), cm[0]) + get_rel(result[2].back(), cm[1]);
            //	return std::abs(result[1].back() - cm[0]) + std::abs(result[2].back()- cm[1]);
        }
        return get_rel(result[3].back(), cm[0]) + get_rel(result[4].back(), cm[1]);
        //	return std::abs(result[3].back() - cm[0]) + std::abs(result[4].back() - cm[1]);

    };

    auto orbit = [&](const bool is_x) {
        notk::NOTKController<double, double> NOTKController;

        NOTKController.add_problem_config(notk_config);

        if (is_x)
        {
            NOTKController.set_borders({-x_max, -1.0}, {x_max, 1.0});
        }
        else
        {
            NOTKController.set_borders({-y_max, -1.0}, {y_max, 1.0});
        }
        NOTKController.set_x_0({0, 0});

        std::vector<std::string> status;

        NOTKController.set_fitness(std::bind(fitness, std::placeholders::_1, is_x));

        if (!NOTKController.check_configuration())
        {
            errorMsg = "An error occured during notk library initialization. Please "
                       "see log files for additional info.";
            return;
        }

        bool flagAbort = true;

        auto result = NOTKController.process(flagAbort, status);

        if (!result)
        {
            errorMsg = "An error occured during notk library calcualtion. Please "
                       "see log files for additional info.";
            return;
        }
        if (!is_x)
        {
            ret_result.first.first = result->get_last_argument();
        }
        else
        {
            ret_result.first.second = result->get_last_argument();
        }
    };
    orbit(false);
    orbit(true);*/

    return ret_result;
}