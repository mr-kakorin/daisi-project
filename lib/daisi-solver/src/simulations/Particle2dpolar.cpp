#include "Particle.h"
#include "Dmath.h"
#include <Constants.h>


template <class PointType>
template <class Archive>
void Particles2dpolar<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
}

template <class PointType>
template <class Archive>
void Particles2dpolar<PointType>::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);

    if (this->positions[0].size() != 0)
        Dmath::Polar2Cartesian(&this->positions[0][0], &this->positions[1][0],
                               &this->positions[2][0], &this->positions[3][0],
                               this->positions[0].size());
}

template <class PointType>
void Particles2dpolar<PointType>::SetCartesianMomentum(int i, const std::vector<PointType>& PX)
{
    this->momentums[0][i] = PX[0] * std::cos(2 * PI() - this->positions[1][i]) +
                            PX[1] * std::sin(2 * PI() - this->positions[1][i]);
    this->momentums[1][i] = PX[0] * std::sin(2 * PI() - this->positions[1][i]) -
                            PX[1] * std::cos(2 * PI() - this->positions[1][i]);
    this->momentums[1][i] = this->momentums[1][i] * this->positions[0][i];
};
template <class PointType>
PointType Particles2dpolar<PointType>::GetCartesianPX(int i)
{
    return this->momentums[0][i] * std::cos(2 * PI() - this->positions[1][i]) +
           (this->momentums[1][i] / this->positions[0][i]) *
               std::sin(2 * PI() - this->positions[1][i]);
};
template <class PointType>
PointType Particles2dpolar<PointType>::GetCartesianPY(int i)
{
    return this->momentums[0][i] * std::sin(2 * PI() - this->positions[1][i]) -
           (this->momentums[1][i] / this->positions[0][i]) *
               std::cos(2 * PI() - this->positions[1][i]);
};

template <class PointType>
void Particles2dpolar<PointType>::SetCartesianPosition(int i, const std::vector<PointType>& xin)
{
    Dmath::Cartesian2Polar(xin[0], xin[1], this->positions[0][i], this->positions[1][i]);
};
template <class PointType>
PointType* Particles2dpolar<PointType>::GetPointerToCartesianX()
{
    return &this->positions[2][0];
};
template <class PointType>
PointType* Particles2dpolar<PointType>::GetPointerToCartesianY()
{
    return &this->positions[3][0];
};

template <class PointType>
void Particles2dpolar<PointType>::PeriodicalEvent()
{
    double dphi    = 0.058177641733144;
    double border1 = 1.5 * PI() - dphi / 2;
    double border2 = 1.5 * PI() + dphi / 2;

    std::vector<PointType> tmp = this->positions[1];

    this->positions[1] = this->positions[4];

    for (int i = 0; i < this->q.size(); i++)
    {

        this->Get_isPeriodical()[i] = 0;

        while (this->positions[1][i] < border1)
            this->positions[1][i] = this->positions[1][i] + dphi;

        while (this->positions[1][i] > border2)
            this->positions[1][i] = this->positions[1][i] - dphi;

        if (this->positions[1][i] < 1.5 * PI())
            this->positions[1][i] =
                1.5 * PI() + (1.5 * PI() - this->positions[1][i]);

        if (!(this->positions[1][i] > 1.5 * PI() && this->positions[1][i] < border2))
        {
            int tt = 0;
        };
    }

    for (int i = 0; i < this->q.size(); i++)
    {
        if ((tmp[i] > this->positions[1][i]) && this->momentums[1][i] > 0)
        {
            this->Get_isPeriodical()[i] = 2;
        }
        if ((tmp[i] < this->positions[1][i]) && this->momentums[1][i] < 0)
        {
            this->Get_isPeriodical()[i] = 1;
        }
    }
};

template <class PointType>
PointType Particles2dpolar<PointType>::GetBeta(int number)
{
    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           (this->momentums[1][number] / this->positions[0][number]) *
                               (this->momentums[1][number] / this->positions[0][number]));
    PointType beta = sqrt(gamma * gamma - 1) / gamma;
    return beta;
}

template <class PointType>
void Particles2dpolar<PointType>::GammaCalc(std::vector<PointType>& gamma)
{
    gamma.resize(this->momentums[0].size());

    for (int number = 0; number < this->momentums[0].size(); number++)
    {

        gamma[number] = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                             (this->momentums[1][number] / this->positions[0][number]) *
                                 (this->momentums[1][number] / this->positions[0][number]));
    }
};

template <class PointType>
void Particles2dpolar<PointType>::GetBetaComponents(PointType& beta1, PointType& beta2,
                                                    PointType& beta3, int number)
{

    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           (this->momentums[1][number] / this->positions[0][number]) *
                               (this->momentums[1][number] / this->positions[0][number]));
    beta1 = this->momentums[0][number] / gamma;
    beta2 = this->momentums[1][number] / (gamma * this->positions[0][number]);
    beta3 = 0;
};

template <class PointType>
PointType Particles2dpolar<PointType>::GetEnergy(int number, PointType mass)
{
    PointType en = LIGHT_VELOCITY() * LIGHT_VELOCITY() * mass;

    PointType gamma = sqrt(1 + this->momentums[0][number] * this->momentums[0][number] +
                           (this->momentums[1][number] / this->positions[0][number]) *
                               (this->momentums[1][number] / this->positions[0][number]));
    return (-(gamma - 1) * en);
};

template <class PointType>
Particles2dpolar<PointType>::Particles2dpolar(int currentSolverType)
{
    this->momentums.resize(2);
    this->coloumbFields.resize(2);
    this->positions.resize(7);
    this->additionalInf.resize(1);
    if (currentSolverType == 1)
    {
        this->additionalInfType.resize(1);
        this->additionalInf.resize(2);
    }
};
template <class PointType>
void Particles2dpolar<PointType>::GetParticlesCloud(int flag, std::vector<void*>& pointArray,
                                                    int& sizeArray, int& sizeElement)
{
    if (this->q.size() <= 0)
    {
        sizeArray   = 0;
        sizeElement = 8;
        return;
    }
    pointArray.resize(2);
    sizeArray     = this->q.size();
    pointArray[0] = (void*)(&this->positions[5][0]);
    pointArray[1] = (void*)(&this->positions[6][0]);
    sizeElement   = sizeof(this->positions[5][0]);
};

template <class PointType>
void Particles2dpolar<PointType>::GetBeamMeasuriments(std::vector<std::vector<PointType>>& data,
                                                      double mass, double Zav, double betaAv){

};

template <class PointType>
void Particles2dpolar<PointType>::GetEmittanceData(std::vector<std::vector<float>>& data,
                                                   int emFlag, double mass, double lambda){};

//// 0 r 1 phi 2 x 3 y 4 phiReal 5 xReal 6 yReal,
////
template class Particles2dpolar<double>;
template class Particles2dpolar<float>;

template void
Particles2dpolar<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);
template void Particles2dpolar<double>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void Particles2dpolar<float>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void
Particles2dpolar<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);