/*#define _USE_MATH_DEFINES
#define SEED    1
#define BRNG    VSL_BRNG_MCG31*/
#include "../base/AccelFlow.h"
#include "FlagStringsSolver.h"
#include "Tools.h"

#include <common_tools/constants.h>

#include <random>
std::vector<std::string> namesFlowParametersSynchrotron = {"Charge number",
                                                           "Mass number",
                                                           "Number of particles",
                                                           "Average Energy, eV",
                                                           "Momentum spread",
                                                           "EmittanceX (norm), pi*mm*mrad",
                                                           "X max, m",
                                                           "dX/dZ max, rad",
                                                           "EmittanceY (norm), pi*mm*mrad",
                                                           "Y max, m",
                                                           "dY/dZ max, rad",
                                                           "Center mass pos, X, m",
                                                           "Center mass pos, Y, m",
                                                           "Center mass pos, dX/dZ, rad",
                                                           "Center mass pos, dY/dZ, rad",
                                                           "Impulse current, A",
                                                           "Channel rel. radius",
                                                           "Distribution type"};
std::vector<std::string> namesFlowParametersRFQ = {"Voltage",
                                                   "Charge number",
                                                   "Mass number",
                                                   "Number of particles",
                                                   "Average Energy, eV",
                                                   "Momentum spread",
                                                   "EmittanceX (norm), pi*mm*mrad",
                                                   "X max, m",
                                                   "dX/dZ max, rad",
                                                   "EmittanceY (norm), pi*mm*mrad",
                                                   "Y max, m",
                                                   "dY/dZ max, rad",
                                                   "Impulse current, A",
                                                   "Channel rel. radius",
                                                   "Distribution type"};

std::vector<std::string> namesFlowParametersDTL = {"Voltage",
                                                   "Charge number",
                                                   "Mass number",
                                                   "Number of particles",
                                                   "Average Energy, eV",
                                                   "Momentum spread",
                                                   "Average phase, rad",
                                                   "Phase spread, rad",
                                                   "EmittanceX (norm), pi*mm*mrad",
                                                   "X max, m",
                                                   "dX/dZ max, rad",
                                                   "EmittanceY (norm), pi*mm*mrad",
                                                   "Y max, m",
                                                   "dY/dZ max, rad",
                                                   "Impulse current, A",
                                                   "Channel rel. radius",
                                                   "Distribution type"};

std::default_random_engine             generator;
std::normal_distribution<double>       distribution(0, 1);
std::uniform_real_distribution<double> distributionUn(0, 1);

template void
AccelFlowBase::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                          const unsigned int file_version);
template void
AccelFlowBase::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                          const unsigned int file_version);

template void
AccelFlowBase::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                     const unsigned int               file_version);

template void
AccelFlowBase::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                     const unsigned int file_version) const;

template <class Archive>
void AccelFlowBase::save(Archive& ar, const unsigned int) const
{
    ar& allParameters;
}
template <class Archive>
void AccelFlowBase::load(Archive& ar, const unsigned int)
{
    ar& allParameters;

    translateParameters();
}

template void
AccelFlow::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                      const unsigned int file_version);
template void
AccelFlow::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                      const unsigned int file_version);

template void AccelFlow::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);

template void
AccelFlow::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                 const unsigned int file_version) const;

template <class Archive>
void AccelFlow::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelFlowBase>(*this);
}
template <class Archive>
void AccelFlow::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelFlowBase>(*this);
}

template void
SynchrotronFlow::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                            const unsigned int file_version);
template void
SynchrotronFlow::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                            const unsigned int file_version);

template void
SynchrotronFlow::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                       const unsigned int file_version);

template void
SynchrotronFlow::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                       const unsigned int file_version) const;

template <class Archive>
void SynchrotronFlow::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelFlowBase>(*this);
}
template <class Archive>
void SynchrotronFlow::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelFlowBase>(*this);
    GenerateParticlesAccel();
}

double LinacDynamicsAccel::CalculateTransmission(double Xaperture, double Yaperture,
                                                 double ChannelRel)
{
    double tr = 0;
    for (int i = 0; i < x.size(); i++)
    {
        if (!isInPerture[i])
            continue;
        if ((x[i] * x[i]) / (Xaperture * Xaperture) + (y[i] * y[i]) / (Yaperture * Yaperture) <
            ChannelRel)
            tr++;
        else
            isInPerture[i] = 0;
    };
    return tr / double(x.size());
};

double SynchrotronDynamics::CalculateTransmission(double Xaperture, double Yaperture,
                                                  double ChannelRel)
{
    double tr = 0;
    for (int i = 0; i < position.size(); i++)
    {
        if (!isInPerture[i])
            continue;
        if ((position[i](0) * position[i](0)) / (Xaperture * Xaperture) +
                (position[i](2) * position[i](2)) / (Yaperture * Yaperture) <
            ChannelRel)
            tr++;
        else
            isInPerture[i] = 0;
    };
    return tr / double(position.size());
};

void SynchrotronFlow::GetMassCenterVector(std::vector<double>& result)
{
    result.resize(5);
    result[0] = centerMassX;
    result[1] = centerMassdX;
    result[2] = centerMassY;
    result[3] = centerMassdY;
    result[4] = momentumSpread;
}

void SynchrotronFlow::SetMassCenterVector(const std::vector<double>& result, const bool is_x)
{
	if (is_x)
	{
		allParameters->set("Center mass pos, X, m", result[0]);
		allParameters->set("Center mass pos, dX/dZ, rad", result[1]);
	}
	else
	{
		allParameters->set("Center mass pos, Y, m", result[0]);
		allParameters->set("Center mass pos, dY/dZ, rad", result[1]);
	}
	translateParameters();
}
void AccelFlowBase::GetParametersAccelFlow(std::vector<std::string>& keysO, std::vector<double>& p)
{
    keysO = allParameters->GetKeys();
    p     = allParameters->GetValues();
};
void AccelFlowBase::SetParametersAccelFlow(const std::vector<double>& in)
{
    allParameters->SetValues(in);
    translateParameters();
};
AccelFlowBase::AccelFlowBase(int ProblemType)
{
    allParameters = new myunsorted_map();
    switch (ProblemType)
    {
    case 0: // synchrotron
        for (int i = 0; i < namesFlowParametersSynchrotron.size(); i++)
            allParameters->insert(namesFlowParametersSynchrotron[i], 0.0);

        break;
    case 1: // rfq
        for (int i = 0; i < namesFlowParametersRFQ.size(); i++)
            allParameters->insert(namesFlowParametersRFQ[i], 0.0);
        break;

    case 3: // dtl
        for (int i = 0; i < namesFlowParametersDTL.size(); i++)
            allParameters->insert(namesFlowParametersDTL[i], 0.0);

        break;
    }
    translateParameters();
};

//#include "mkl_vsl.h"
//
/*
void LinacFlow::GenerateOneParticleRFQ(double lambda, double z0)
{
        oneParticle.resize(8);
        double restEnergy = -mass*commtools::LIGHT_VELOCITY()*commtools::LIGHT_VELOCITY() /
commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

        double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
        double beta = sqrt(gamma * gamma - 1) / gamma;

        double P = gamma*beta;
        oneParticle[0] = P;
        oneParticle[1] = z0;

        double emittanceXL = emittanceX*1e-6;
        double emittanceYL = emittanceY*1e-6;

        double dXmaxL = dXmax*beta;

        double betaX = Xmax*Xmax / emittanceXL;
        double gammaX = dXmaxL*dXmaxL / emittanceXL;
        double betaY = Xmax*Xmax / emittanceYL;
        double gammaY = dXmaxL*dXmaxL / emittanceYL;
        double AlphaX = sqrt(gammaX*betaX - 1);
        if (dXmax < 0)
                AlphaX = -AlphaX;

        if (std::abs(gammaX*betaX - 1) < 1e-7)
                AlphaX = 0;

        double AlphaY = sqrt(gammaY*betaY - 1);
        if (dXmax < 0)
                AlphaY = -AlphaY;

        if (std::abs(gammaY*betaY - 1) < 1e-7)
                AlphaY = 0;

        oneParticle[2] = emittanceXL*betaX;
        oneParticle[3] = -emittanceXL*AlphaX;
        oneParticle[4] = emittanceXL*gammaX;

        oneParticle[5] = emittanceYL*betaY;
        oneParticle[6] = -emittanceYL*AlphaY;
        oneParticle[7] = emittanceYL*gammaY;

        oneParticleData.resize(4);
        oneParticleData[0] = -2;


}
*/

std::string AccelFlowBase::getParticleType()
{
    if (massNumber == 1 && chargeNumber == 1)
        return "PROTON";

    if (massNumber == 197)
        return "AURUM";
};
double AccelFlowBase::getRestMassInGeV()
{
    return mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
           std::abs(commtools::ELECTRON_CHARGE()) * 1e-9;
};

std::vector<void*> LinacDynamicsLong::GetData()
{
    std::vector<void*> result(2);

    if (p.size() == 0)
        return result;

    result[0] = (void*)(&z[0]);
    result[1] = (void*)(&p[0]);

    return result;
}

std::vector<void*> LinacDynamicsAccel::GetData()
{
    std::vector<void*> result(6);

    if (z.size() == 0)
        return result;

    result[0] = (void*)(&x[0]);
    result[1] = (void*)(&dx[0]);
    result[2] = (void*)(&y[0]);
    result[3] = (void*)(&dy[0]);
    result[4] = (void*)(&z[0]);
    result[5] = (void*)(&dz[0]);

    return result;
}

std::vector<void*> LinacDynamicsTransv::GetData()
{
    std::vector<void*> result(8);

    if (p.size() == 0)
        return result;

    result[0] = (void*)(&z[0]);
    result[1] = (void*)(&p[0]);
    result[2] = (void*)(&S11x[0]);
    result[3] = (void*)(&S12x[0]);
    result[4] = (void*)(&S22x[0]);
    result[5] = (void*)(&S11y[0]);
    result[6] = (void*)(&S12y[0]);
    result[7] = (void*)(&S22y[0]);

    return result;
}

double AccelFlow::GetOutVelocity() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(outEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta;
}
/*
void AccelFlow::GenerateParticlesRFQ(double lambda, double z0, double Rchann)
{
        dynamicsLong.clearandResize(nParticlesLong);
        dynamicsTransv.clearandResize(nParticlesTransv);

        double restEnergy = -mass*commtools::LIGHT_VELOCITY()*commtools::LIGHT_VELOCITY() /
commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

        double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
        double beta = sqrt(gamma * gamma - 1) / gamma;

        double P = gamma*beta;

        VSLStreamStatePtr stream;

        vslNewStream(&stream, BRNG, SEED);

        vdRngGaussian(0, stream, nParticlesLong-1, &dynamicsLong.p[1], P, P*momentumSpread / 3);

        vdRngUniform(0, stream, nParticlesLong-1, &dynamicsLong.z[1], -lambda*beta + z0, z0);

        dynamicsLong.p[0] = P;
        dynamicsLong.z[0] = -lambda*beta/2 + z0;

        double R = channelRelRadius*Rchann;
        dynamicsLong.longQ = impulseCurrent*lambda /
(commtools::LIGHT_VELOCITY()*nParticlesLong*commtools::PI()*R*R);

        if (nParticlesTransv == 1)
                dynamicsLong.z[0] = z0;

        for (int i = 0; i < nParticlesLong; i++)
        {
                dynamicsLong.cellNumber[i] = -2;
        };



        double emittanceXL = emittanceX*1e-6;
        double emittanceYL = emittanceY*1e-6;

        double dXmaxL = dXmax*beta;

        double betaX = Xmax*Xmax / emittanceXL;
        double gammaX = dXmaxL*dXmaxL / emittanceXL;
        double betaY = Xmax*Xmax / emittanceYL;
        double gammaY = dXmaxL*dXmaxL / emittanceYL;
        double AlphaX = sqrt(gammaX*betaX - 1);
        if (dXmax < 0)
                AlphaX = -AlphaX;

        if (std::abs(gammaX*betaX - 1) < 1e-7)
                AlphaX = 0;

        double AlphaY = sqrt(gammaY*betaY - 1);
        if (dXmax < 0)
                AlphaY = -AlphaY;

        if (std::abs(gammaY*betaY - 1) < 1e-7)
                AlphaY = 0;

        double dZ = lambda*beta / nParticlesTransv;

        for (int i = 0; i < nParticlesTransv; i++)
        {
                dynamicsTransv.S11x[i] = emittanceXL*betaX;
                dynamicsTransv.S12x[i] = -emittanceXL*AlphaX;
                dynamicsTransv.S22x[i] = emittanceXL*gammaX;

                dynamicsTransv.S11y[i] = emittanceYL*betaY;
                dynamicsTransv.S12y[i] = -emittanceYL*AlphaY;
                dynamicsTransv.S22y[i] = emittanceYL*gammaY;
                dynamicsTransv.cellNumber[i] = -2;
                dynamicsTransv.z[i] = -lambda*beta + i*dZ + dZ / 2 + z0;
        };

        vdRngGaussian(0, stream, nParticlesTransv - 1, &dynamicsTransv.p[1], P, P*momentumSpread /
3);

        //vdRngUniform(0, stream, nParticlesTransv - 1, &dynamicsTransv.z[1], -lambda*beta + z0,
z0);
        //dynamicsTransv.p[i] = P;
//	dynamicsTransv.z[i] = -lambda*beta + i*dZ + dZ / 2 + z0;

        dynamicsTransv.p[0] = P;
        dynamicsTransv.z[0] = -lambda*beta / 2 + z0;

        dynamicsTransv.transvQ = impulseCurrent*lambda /
(commtools::LIGHT_VELOCITY()*nParticlesTransv);

        if (nParticlesTransv == 1)
                dynamicsTransv.z[0] = z0;
};*/

double AccelFlow::GetStartVelocity() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta;
};
double AccelFlow::GetStartMomentum() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta * gamma;
}
double AccelFlow::GetAlpha() const
{
    return charge / (mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY());
};
double AccelFlow::GetRestEnergy() const
{
    return mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY();
};

void AccelFlowBase::translateParameters()
{
    allParameters->find("Voltage", voltage);

    allParameters->find("Mass number", massNumber);

    allParameters->find("Charge number", chargeNumber);

    allParameters->find("EmittanceX (norm), pi*mm*mrad", emittanceX);

    allParameters->find("EmittanceY (norm), pi*mm*mrad", emittanceY);

    allParameters->find("EmittanceX (unnorm), pi*mm*mrad", emittanceXunnorm);

    allParameters->find("EmittanceY (unnorm), pi*mm*mrad", emittanceYunnorm);
    phaseSpread = 0;
    allParameters->find("Average phase, rad", averagePhase);
    allParameters->find("Phase spread, rad", phaseSpread);

    double averagePhase;
    double phaseSpread;

    allParameters->find("Average Energy, eV", averageEnergy);

    allParameters->find("Momentum spread", momentumSpread);

    allParameters->find("Impulse current, A", impulseCurrent);

    double tmp;
    allParameters->find("Number of particles", tmp);
    nParticles = tmp;
    allParameters->find("Number of particles long", tmp);
    nParticlesLong = tmp;
    allParameters->find("Number of particles transv", tmp);
    nParticlesTransv = tmp;

    allParameters->find("X max, m", Xmax);

    allParameters->find("dX/dZ max, rad", dXmax);

    allParameters->find("Y max, m", Ymax);

    allParameters->find("dY/dZ max, rad", dYmax);

    allParameters->find("Out energy", outEnergy);

    allParameters->find("Channel rel. radius", channelRelRadius);
    allParameters->find("Center mass pos, X, m", centerMassX);
    allParameters->find("Center mass pos, Y, m", centerMassY);
    allParameters->find("Center mass pos, dX/dZ, rad", centerMassdX);
    allParameters->find("Center mass pos, dY/dZ, rad", centerMassdY);
    allParameters->find("Distribution type", tmp);
    distrType = tmp;

    charge = -chargeNumber * commtools::ELECTRON_CHARGE();
    mass   = massNumber * commtools::PROTON_MASS();
};

void LinacDynamicsLong::clearandResize(int nParticles)
{
    NparticlesLong = nParticles;
    p.clear();
    z.clear();

    Ez.clear();
    cellNumber.clear();

    p.resize(nParticles);
    z.resize(nParticles);
    Ez.resize(nParticles);
    cellNumber.resize(nParticles);
};

void LinacDynamicsTransv::clearandResize(int nParticles)
{
    NparticlesTransv = nParticles;
    p.clear();
    z.clear();
    S11x.clear();
    S12x.clear();
    S22x.clear();
    S11y.clear();
    S12y.clear();
    S22y.clear();
    Ex.clear();
    Ey.clear();
    Ez.clear();
    cellNumber.clear();

    p.resize(nParticles);
    z.resize(nParticles);
    S11x.resize(nParticles);
    S12x.resize(nParticles);
    S22x.resize(nParticles);
    S11y.resize(nParticles);
    S12y.resize(nParticles);
    S22y.resize(nParticles);
    Ex.resize(nParticles);
    Ey.resize(nParticles);
    Ez.resize(nParticles);
    EzCol.resize(nParticles);
    cellNumber.resize(nParticles);
};

void LinacDynamicsAccel::clearandResize(int nParticles)
{
    int    Nparticles;
    double longQ;

    Nparticles = nParticles;
    zStart.clear();
    dzStart.clear();
    xStart.clear();
    dxStart.clear();
    yStart.clear();
    dyStart.clear();
    cellNumber.clear();
    isInPerture.clear();

    zStart.resize(nParticles);
    dzStart.resize(nParticles);
    xStart.resize(nParticles);
    dxStart.resize(nParticles);
    yStart.resize(nParticles);
    dyStart.resize(nParticles);
    cellNumber.resize(nParticles);
    isInPerture.resize(nParticles);
}

void SynchrotronDynamics::clearandResize(int nParticles)
{
    Nparticles = nParticles;

    isInPerture.clear();
    position.clear();
    positionStart.clear();

    isInPerture.resize(nParticles);
    position.resize(nParticles);
    positionStart.resize(nParticles);
}

void LinacDynamicsAccel::searchBadParticle(std::vector<unsigned int>& Indexes)
{
    for (int i = 0; i < isInPerture.size(); i++)
    {
        if (0 == isInPerture[i])
        {
            int flagAdd = 0;
            for (int j = 0; j < Indexes.size(); j++)
            {
                if (Indexes[j] == i)
                {
                    flagAdd = 1;
                    break;
                };
            }
            if (0 == flagAdd)
                Indexes.push_back(i);
        }
    };
};

void LinacDynamicsAccel::removeParticle()
{
    //	searchBadParticle(indexes);

    int n1 = int(zStart.size());
    int j;
    for (int i = 0; i < n1; i++)
    {
        if (!isInPerture[i])
        {
            j = i;
            while (!isInPerture[j])
            {
                j++;
                if (j >= n1)
                    break;
            }

            zStart.erase(zStart.begin() + i, zStart.begin() + j);
            dzStart.erase(dzStart.begin() + i, dzStart.begin() + j);
            xStart.erase(xStart.begin() + i, xStart.begin() + j);
            dxStart.erase(dxStart.begin() + i, dxStart.begin() + j);
            yStart.erase(yStart.begin() + i, yStart.begin() + j);
            dyStart.erase(dyStart.begin() + i, dyStart.begin() + j);
            isInPerture.erase(isInPerture.begin() + i, isInPerture.begin() + j);

            n1 = n1 - (j - i);
            i  = i - 1;
        };
    }
};

void SynchrotronDynamics::removeParticle()
{
    //	searchBadParticle(indexes);

    int n1 = int(positionStart.size());
    int j;
    for (int i = 0; i < n1; i++)
    {
        if (!isInPerture[i])
        {
            j = i;
            while (!isInPerture[j])
            {
                j++;
                if (j >= n1)
                    break;
            }

            positionStart.erase(positionStart.begin() + i, positionStart.begin() + j);
            isInPerture.erase(isInPerture.begin() + i, isInPerture.begin() + j);

            n1 = n1 - (j - i);
            i  = i - 1;
        };
    }
};

void LinacDynamicsAccel::removeParticle(std::vector<unsigned int>& indexes)
{
    //	searchBadParticle(indexes);

    int n      = int(indexes.size());
    Nparticles = Nparticles - n;
    for (int i = 0; i < n; i++)
    {
        zStart[indexes[i]] = -1;
    };

    int n1 = int(zStart.size());
    int j;
    for (int i = 0; i < n1; i++)
    {
        if (zStart[i] == -1)
        {
            j = i;
            while (zStart[j] == -1)
            {
                j++;
                if (j >= n1)
                    break;
            }

            zStart.erase(zStart.begin() + i, zStart.begin() + j);
            dzStart.erase(dzStart.begin() + i, dzStart.begin() + j);
            xStart.erase(xStart.begin() + i, xStart.begin() + j);
            dxStart.erase(dxStart.begin() + i, dxStart.begin() + j);
            yStart.erase(yStart.begin() + i, yStart.begin() + j);
            dyStart.erase(dyStart.begin() + i, dyStart.begin() + j);
            isInPerture.erase(isInPerture.begin() + i, isInPerture.begin() + j);

            z.erase(z.begin() + i, z.begin() + j);
            dz.erase(dz.begin() + i, dz.begin() + j);
            x.erase(x.begin() + i, x.begin() + j);
            dx.erase(dx.begin() + i, dx.begin() + j);
            y.erase(y.begin() + i, y.begin() + j);
            dy.erase(dy.begin() + i, dy.begin() + j);

            n1 = n1 - (j - i);
            i  = i - 1;
        };
    }
};

int LinacDynamicsTransv::checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                                           const std::vector<double>& LL,
                                           const std::vector<double>& MinumalRadii, double lambda,
                                           double relChann)
{
    double zAv = 0;
    double pAv = 0;

    int    n    = 0;
    double maxR = 0;
    double Lend = LL.back();
    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    for (int i = 0; i < z.size(); i++)
    {
        if (cellNumber[i] != -1)
            n++;
    };
    for (int i = 0; i < z.size(); i++)
    {
        //	if (cellNumber[i] != -1)
        //	{
        double R = std::max(S11x[i], S11y[i]);
        if (R > maxR)
            maxR = R;
        zAv      = zAv + z[i];
        pAv      = pAv + p[i];
        //	};
    };
    zAv = zAv / z.size();
    pAv = pAv / z.size();

    cell = 0;
    for (cell = 0; cell < LL.size(); cell++)
    {
        if (LL[cell] > z[0])
        {
            cell--;
            break;
        }
    };
    if (cell == -1)
        cell = 0;
    if (cell >= LL.size() - 1)
        cell = LL.size() - 2;

    double Rchannel = MinumalRadii[cell];
    OutParams[0].push_back(Rchannel);
    OutParams[1].push_back(sqrt(maxR));

    double transmission = 0;
    double acceleration = 0;

    for (int i = 0; i < z.size(); i++)
    {
        // if (z[i]>Lend || z[i]<0)
        //	cellNumber[i] = -1;

        if (z[i] > Lend)
        {
            cellNumber[i] = -1;
        }
        double R = std::max(sqrt(S11x[i]), sqrt(S11y[i]));
        if (S11x[i] < 1e-8 || S11y[i] < 1e-8)
            R = 1.0000e-04;

		if (R > 0.0058)
		{
			int tt = 0;
		}

        if (R < Rchannel * relChann || cell < 20)
            transmission++;
        else
            cellNumber[i] = -1;

        double dPh = 2 * commtools::PI() * (zAv - z[i]) / (pAv * lambda);

        if (dPh < commtools::PI())
            acceleration++;
    }

    if (acceleration == 0)
    {
        int tt = 0;
    }
    OutParams[3].push_back(acceleration / z.size());

    OutParams[2].push_back(transmission / z.size());
    OutParams[4].push_back(pAv);

    return n;
};

int LinacDynamicsLong::checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                                         const std::vector<double>& LL,
                                         const std::vector<double>& MinumalRadii, double lambda)
{
    double zAv = z[0];
    double pAv = 0;

    int    n    = 0;
    double maxR = 0;
    double Lend = LL.back();
    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    for (int i = 0; i < z.size(); i++)
    {
        if (cellNumber[i] != -1)
            n++;
    };
    for (int i = 0; i < z.size(); i++)
    {
        //		zAv = zAv + z[i];
        pAv = pAv + p[i];
    };
    //	zAv = zAv / z.size();
    pAv = pAv / z.size();

    cell = 0;
    for (cell = 0; cell < LL.size(); cell++)
    {
        if (LL[cell] > zAv)
        {
            cell--;
            break;
        }
    };
    if (cell == -1)
        cell = 0;
    if (cell >= LL.size() - 1)
        cell = LL.size() - 2;

    double transmission = 0;
    double acceleration = 0;

    for (int i = 0; i < z.size(); i++)
    {
        //	if (z[i]>Lend || z[i]<0)
        //		cellNumber[i] = -1;

        if (z[i] > Lend)
            cellNumber[i] = -1;
        double dPh        = 2 * commtools::PI() * (zAv - z[i]) / (pAv * lambda);

        if (dPh < commtools::PI())
            acceleration++;

        //	if (std::abs(pAv - p[i]) / pAv < 0.04)
        //	acceleration++;
    }

    if (acceleration == 0)
    {
        int tt = 0;
    }

    OutParams[3].push_back(acceleration / z.size());

    OutParams[4].push_back(pAv);

    return n;
};

void SynchrotronFlow::GenerateParticlesAccelTest(int nparticles, double xMax, double yMax,
                                                 double dxMax, double dyMax)
{
    dynamicsAccelTest.clearandResize(nparticles);

    for (int i = 0; i < nparticles; i++)
    {
        dynamicsAccelTest.positionStart[i].resize(5);

        dynamicsAccelTest.positionStart[i](0) = -xMax + distributionUn(generator) * (2 * xMax);
        dynamicsAccelTest.positionStart[i](2) = -yMax + distributionUn(generator) * (2 * yMax);
        dynamicsAccelTest.positionStart[i](1) = -dxMax + distributionUn(generator) * (2 * dxMax);
        dynamicsAccelTest.positionStart[i](3) = -dyMax + distributionUn(generator) * (2 * dyMax);
        dynamicsAccelTest.positionStart[i](4) = 0;
    }
};
void AccelFlowBase::getTransvData(std::vector<double>& x, std::vector<double>& dx,
                                  std::vector<double>& y, std::vector<double>& dy)
{

    x.resize(nParticles);
    dx.resize(nParticles);
    y.resize(nParticles);
    dy.resize(nParticles);

    std::default_random_engine             generator;
    std::normal_distribution<double>       distribution(0, 1);
    std::uniform_real_distribution<double> distributionUn(0, 1);

    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    double P = gamma * beta;

    emittanceXunnorm = emittanceX / beta;
    emittanceYunnorm = emittanceY / beta;

    double Xparam[3];
    double Yparam[3];
    double Angles[3];
    double e[2];

    double alpha[3];

    Xparam[0] = Xmax;                    // double X
    Yparam[0] = dXmax;                   // dX
    e[0]      = emittanceXunnorm * 1e-6; // em

    Xparam[1] = Ymax;                    // double Y
    Yparam[1] = dYmax;                   // dY
    e[1]      = emittanceYunnorm * 1e-6; // em

    double R = std::max(Xparam[2], Xparam[1]);

    double gammaE[2];
    double betaE[2];
    double phi[2];
    double k[2];

    double Major[2];
    double Minor[2];

    for (int i = 0; i < 2; i++)
    {
        //	e[i] = Xparam[i] * Yparam[i] / sqrt(1 + alpha[i] * alpha[i]);

        alpha[i]  = sqrt(Xparam[i] * Xparam[i] * Yparam[i] * Yparam[i] - e[i] * e[i]) / e[i];
        gammaE[i] = Yparam[i] * Yparam[i] / std::abs(e[i]);
        betaE[i]  = Xparam[i] * Xparam[i] / std::abs(e[i]);

        phi[i] = atan(2 * alpha[i] / (gammaE[i] - betaE[i])) / 2;

        k[i] = tan(phi[i]);

        double x = sqrt(std::abs(e[i]) / (gammaE[i] + 2 * alpha[i] * k[i] + betaE[i] * k[i] * k[i]));
        double y = k[i] * x;
        Major[i] = sqrt(x * x + y * y);
        Minor[i] = std::abs(e[i]) / (Major[i]);

        if (Xparam[i] > 0)
            phi[i] = phi[i] + commtools::PI() / 2;
    };

    // double EnergyDistr(particlesPerBunch);
    // double PhaseDistr[particlesPerBunch];
    // double r_Distr[particlesPerBunch];
    // double pr_Distr[particlesPerBunch];
    std::vector<std::vector<double>> Data(6);

    auto generate = [&](const double Major, const double Minor, const double phi,
                        const double centerMassX,
                        const double centerMassdX) -> std::pair<double, double> {
        double xtmp, ytmp;
        std::pair<double, double> result;
        if (distrType == 0)
        {
            xtmp = (Major / 3) * distribution(generator);
            ytmp = (Minor / 3) * distribution(generator);
        };
        if (distrType == 1)
        {
            while (1)
            {
                xtmp = -Major + distributionUn(generator) * (2 * Major);
                ytmp = -Minor + distributionUn(generator) * (2 * Minor);
                if ((xtmp * xtmp / (Major * Major) + ytmp * ytmp / (Minor * Minor)) < 1)
                    break;
            }
        };

        if (Xparam[0] < 0)
        {
            result.first  = xtmp * std::cos(phi) - ytmp * std::sin(phi) + centerMassX;
            result.second = xtmp * std::sin(phi) + ytmp * std::cos(phi) + centerMassdX;
        }
        else
        {
            result.first  = xtmp * std::cos(phi) - ytmp * std::sin(phi) + centerMassdX;
            result.second = xtmp * std::sin(phi) + ytmp * std::cos(phi) + centerMassX;
        }
        return result;
    };

    for (int i = 0; i < nParticles; i++)
    {
        auto xx = generate(Major[0], Minor[0], phi[0], centerMassX, centerMassdX);
        auto yy = generate(Major[1], Minor[1], phi[1], centerMassY, centerMassdY);

        x[i]  = xx.second;
        dx[i] = xx.first;

        y[i]  = yy.second;
        dy[i] = yy.first;
    }
}
void SynchrotronFlow::GenerateParticlesAccel()
{
    dynamicsAccel.clearandResize(nParticles);
    std::vector<double> x;
    std::vector<double> dx;
    std::vector<double> y;
    std::vector<double> dy;
    getTransvData(x, dx, y, dy);
    std::fill(dynamicsAccel.isInPerture.begin(), dynamicsAccel.isInPerture.end(), 1);
    for (size_t i = 0; i < nParticles; i++)
    {
        dynamicsAccel.positionStart[i].resize(5);
        dynamicsAccel.positionStart[i](0) = x[i];
        dynamicsAccel.positionStart[i](1) = dx[i];
        dynamicsAccel.positionStart[i](2) = y[i];
        dynamicsAccel.positionStart[i](3) = dy[i];
        dynamicsAccel.positionStart[i](4) = momentumSpread;
    }
}
void AccelFlow::GenerateParticlesAccel()
{
    dynamicsAccel.clearandResize(nParticles);

    getTransvData(dynamicsAccel.xStart, dynamicsAccel.dxStart, dynamicsAccel.yStart,
                  dynamicsAccel.dyStart);
    std::fill(dynamicsAccel.isInPerture.begin(), dynamicsAccel.isInPerture.end(), 1);
};

void LinacDynamicsAccel::Init()
{
    z  = zStart;
    dz = dzStart;
    x  = xStart;
    dx = dxStart;
    y  = yStart;
    dy = dyStart;
}

void SynchrotronDynamics::Init()
{
    position = positionStart;
}

void SynchrotronFlow::SaveToMadFile(const std::string& fileName)
{
    FILE* fid;
    fid = fopen(fileName.c_str(), "w");
    for (int i = 0; i < nParticles; i++)
    {
        fprintf(fid, "ptc_start, x=%lf, px=%lf, y=%lf, py=%lf;\n",
                dynamicsAccel.positionStart[i](0), dynamicsAccel.positionStart[i](1),
                dynamicsAccel.positionStart[i](2), dynamicsAccel.positionStart[i](3));
    };
    fclose(fid);
};
void SynchrotronFlow::SaveBetaToMadFile(const std::string& fileName)
{
    FILE* fid;
    fid = fopen(fileName.c_str(), "w");
    std::vector<double> Xtwiss;
    std::vector<double> Ytwiss;
    arma::vec           res = GetTwissVector();
    fprintf(fid, "initial: BETA0, BETX=%lf, ALFX=%lf, MUX=%lf, BETY=%lf, ALFY=%lf, MUY=%lf;",
            res(0), res(1), 0.0, res(3), res(4), 0.0);
    fclose(fid);
};
/*void AccelFlowBase::GetTwissVectorX(std::vector<double>& Xtwiss)
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma     = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta      = sqrt(gamma * gamma - 1) / gamma;
    double P         = gamma * beta;
    emittanceXunnorm = emittanceX / beta;
    double e         = emittanceXunnorm * 1e-6;

    double alphae = sqrt(Xmax * Xmax * dXmax * dXmax - e * e) / e;
    if (dXmax < 0 || Xmax < 0)
        alphae = -alphae;
    double gammae = dXmax * dXmax / e;
    double betae  = Xmax * Xmax / e;
    Xtwiss.resize(3);
    Xtwiss[0] = betae;
    Xtwiss[1] = alphae;
    Xtwiss[2] = gammae;
};*/
arma::vec AccelFlowBase::GetTwissVector()
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma     = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta      = sqrt(gamma * gamma - 1) / gamma;
    double P         = gamma * beta;
    emittanceYunnorm = emittanceY / beta;
    emittanceXunnorm = emittanceX / beta;

    double ey = emittanceYunnorm * 1e-6;
    double ex = emittanceXunnorm * 1e-6;

    double alphaey = sqrt(Ymax * Ymax * dYmax * dYmax - ey * ey) / ey;
    double gammaey = dYmax * dYmax / ey;
    double betaey  = Ymax * Ymax / ey;
    if (dYmax < 0 || Ymax < 0)
        alphaey = -alphaey;

    double alphaex = sqrt(Xmax * Xmax * dXmax * dXmax - ex * ex) / ex;
    double gammaex = dXmax * dXmax / ex;
    double betaex  = Xmax * Xmax / ex;
    if (dXmax < 0 || Xmax < 0)
        alphaex = -alphaex;

    arma::vec result(6);
    result(0) = betaex;
    result(1) = alphaex;
    result(2) = gammaex;
    result(3) = betaey;
    result(4) = alphaey;
    result(5) = gammaey;
    return result;
};
double AccelFlowBase::GetTotalEnergy() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    return restEnergy + averageEnergy;
};