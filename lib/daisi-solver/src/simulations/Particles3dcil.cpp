#include "Particle.h"
#include "Dmath.h"
#include <common_tools/constants.h>

template class Particles3dcil<double>;
template class Particles3dcil<float>;

template <class PointType>
void Particles3dcil<PointType>::GetBeamMeasuriments(std::vector<std::vector<PointType>>& data,
                                                    double mass, double Zav, double betaAv){
    /*PointType xTmp;
    PointType yTmp;
    PointType pxTmp;
    PointType pyTmp;
    PointType beta2;

    float en = commtools::LIGHT_VELOCITY()* commtools::LIGHT_VELOCITY()*mass / (1e3 *
commtools::ELECTRON_CHARGE());




    for (int i = 0; i < r.size(); i++)
    {
            float energy = (-(gamma[i] - 1)*en);

            float gammaAv = sqrt(1 + betaAv*betaAv);
            float energyAv = (-(gammaAv - 1)*en);

            Dmath::Polar2Cartesian(r[i], phi[i], xTmp, yTmp);

            if (r[i]>0.006 || std::abs((energy - energyAv) / energyAv) > 0.05)
                    continue;

            beta2 = pz[i] / gamma[i];commtools::PI()commtools::PI()
commtools::PI()commtools::PI()
            pxTmp = (pr[i] * std::cos(2 * Dconst::PI - phi[i]) + (pphi[i] / r[i])* std::sin(2 * Dconst::PI -
phi[i])) / (gamma[i]
    * beta2); pyTmp = (pr[i] * std::sin(2 * Dconst::PI - phi[i]) - (pphi[i] / r[i]) * std::cos(2 * Dconst::PI
- phi[i])) /
    (gamma[i] * beta2);


            data[0][0] = data[0][0] + xTmp*xTmp;
            data[0][1] = data[0][1] + pxTmp*pxTmp;
            data[0][2] = data[0][2] + xTmp;
            data[0][3] = data[0][3] + pxTmp;
            data[0][4] = data[0][4] + xTmp*pxTmp;

            data[1][0] = data[1][0] + yTmp*yTmp;
            data[1][1] = data[1][1] + pyTmp*pyTmp;
            data[1][2] = data[1][2] + yTmp;
            data[1][3] = data[1][3] + pyTmp;
            data[1][4] = data[1][4] + yTmp*pyTmp;



            data[2][0] = data[2][0] + (z[i] - Zav) * (z[i] - Zav);
            data[2][1] = data[2][1] + (energy - energyAv) * (energy - energyAv);
            data[2][2] = data[2][2] + (z[i] - Zav);
            data[2][3] = data[2][3] + (energy - energyAv);
            data[2][4] = data[2][4] + (z[i] - Zav) * (energy - energyAv);

            data[3][0] = data[3][0] + beta2;
    }*/

};

template <class PointType>
void Particles3dcil<PointType>::GetEmittanceData(std::vector<std::vector<float>>& data, int emFlag,
                                                 double mass, double lambda){
    /*data.resize(2);
    data[0].resize(r.size());
    data[1].resize(r.size());
    PointType en = commtools::LIGHT_VELOCITY()* commtools::LIGHT_VELOCITY()*mass;

    double enAv = 0;commtools::PI()
    double zAv = 0;
    double omega = 2 * Dconst::PI / lambda;

    float tmp;
    float beta2;

    switch (emFlag)
    {
    case 0:

            for (int i = 0; i < gamma.size(); i++)
            {
                    data[1][i]=-(gamma[i] - 1)*en;commtools::PI()
                    enAv = enAv + (-(gamma[i] - 1)*en);
                    data[0][i] = 180 * z[i] * omega / Dconst::PI;
                    while (data[0][i]>180)
                            data[0][i] = data[0][i] - 360;
                    zAv = zAv + data[0][i];
            };

            if (gamma.size()!=0)
            {
                    zAv = zAv / gamma.size();
                    enAv = enAv / gamma.size();
            }
            for (int i = 0; i < gamma.size(); i++)
            {
                    data[1][i] = 100*(data[1][i] - enAv) / enAv;
                    data[0][i] = data[0][i] - zAv;
            };
            break;
    case 1:
            for (int i = 0; i < gamma.size(); i++)
            {
                    beta2 = pz[i] / gamma[i];commtools::PI()commtools::PI()
                    Dmath::Polar2Cartesian(float(r[i]), float(phi[i]), data[0][i], tmp);
                    data[1][i] = (pr[i] * std::cos(2 * Dconst::PI - phi[i]) + (pphi[i] / r[i])* std::sin(2 *
    Dconst::PI - phi[i]))
    / (gamma[i]*beta2);
            };
            break;
    case 2:
            for (int i = 0; i < gamma.size(); i++)
            {
                    beta2 = pz[i] / gamma[i];commtools::PI()commtools::PI()
                    Dmath::Polar2Cartesian(float(r[i]), float(phi[i]), tmp, data[0][i]);
                    data[1][i] = (pr[i] * std::sin(2 * Dconst::PI - phi[i]) - (pphi[i] / r[i]) * std::cos(2 *
    Dconst::PI -
    phi[i])) / (gamma[i] * beta2);
            };
            break;
    case 3:
            for (int i = 0; i < gamma.size(); i++)
                    Dmath::Polar2Cartesian(float(r[i]), float(phi[i]), data[0][i], data[1][i]);
            break;
    };*/
};

template <class PointType>
PointType Particles3dcil<PointType>::GetBeta(int number)
{
    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           this->momentums[1][number] * this->momentums[1][number] +
                           (this->momentums[2][number] / this->positions[0][number]) *
                               (this->momentums[2][number] / this->positions[0][number]));
    PointType beta = sqrt(gamma * gamma - 1) / gamma;
    return beta;
}

template <class PointType>
void Particles3dcil<PointType>::GetBetaComponents(PointType& beta1, PointType& beta2,
                                                  PointType& beta3, int number)
{
    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           this->momentums[1][number] * this->momentums[1][number] +
                           (this->momentums[2][number] / this->positions[0][number]) *
                               (this->momentums[2][number] / this->positions[0][number]));

    beta1 = this->momentums[0][number] / gamma;
    beta2 = this->momentums[1][number] / gamma;
    beta3 = this->momentums[2][number] / (gamma * this->positions[0][number]);
};

template <class PointType>
void Particles3dcil<PointType>::GammaCalc(std::vector<PointType>& gamma)
{
    gamma.resize(this->momentums[0].size());

    for (int number = 0; number < this->momentums[0].size(); number++)
    {

        gamma[number] = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                             this->momentums[1][number] * this->momentums[1][number] +
                             (this->momentums[2][number] / this->positions[0][number]) *
                                 (this->momentums[2][number] / this->positions[0][number]));
    }
};

template <class PointType>
PointType Particles3dcil<PointType>::GetEnergy(int number, PointType mass)
{
    PointType en    = commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() * mass;
    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           this->momentums[1][number] * this->momentums[1][number] +
                           (this->momentums[2][number] / this->positions[0][number]) *
                               (this->momentums[2][number] / this->positions[0][number]));
    return (-(gamma - 1) * en);
};

template <class PointType>
Particles3dcil<PointType>::Particles3dcil(int currentSolverType)
{
    this->momentums.resize(3);
    this->positions.resize(3);
    this->coloumbFields.resize(2);
    this->magneticFields.resize(2);
    if (currentSolverType == 1)
    {
        this->additionalInfType.resize(1);
        this->additionalInf.resize(1);
    }
};

/*template <class PointType>
void Particles3dcil<PointType>::SetPartCharge(double I, int frequency)
{
        for (int k = 0; k < q.size(); k++)
        {
                q[k] = I / (NParticles*frequency);
        }
};*/