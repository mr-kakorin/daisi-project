#define _USE_MATH_DEFINES
#define SEED 1
#define BRNG VSL_BRNG_MCG31

#include "stdio.h"
#include <random>

#include <common_tools/constants.h>

#include "Dmath.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice2daxsVectorized.h"
#include "Geom.h"
#include "GridData.h"
#include "Particle.h"
#include "ParticleSource.h"
#include "Tools.h"

static std::default_random_engine       generator;
static std::normal_distribution<double> distribution(0, 1);

double andgleDistribution(double t, double t0)
{
    return exp(-(std::sin(t) * std::sin(t)) / (std::sin(t0) * std::sin(t0)));
};
template class EmitterDevice2daxs<float>;
template class EmitterDevice2daxs<double>;

template <class PointType>
void EmitterDevice2daxs<PointType>::PreliminaryGeneration(
    const std::shared_ptr<Particles3dcil<PointType>>& particlesDataZ0, PointType restMass){
    /*particlesData->clear();
    int particlesPerBunch = this->NumbersParams[0];

    particlesData->resize(particlesPerBunch);
    particlesDataZ0->resize(particlesPerBunch);



    std::vector<double> en (particlesPerBunch);
    for (int j = 0; j < particlesPerBunch; j++)
    {
            Dmath::Cartesian2Polar(PointType(Data[2][j]), PointType(Data[4][j]),
    particlesData->Get_r()[j], particlesData->Get_phi()[j]);
    particlesData->SetCartesianMomentumPolar(j, PointType(Data[3][j] * beta), PointType(Data[5][j] *
    beta));

            particlesData->Get_z()[j] = (Data[0][j] + this->DistribParams[4]) * M_PI / 180;

            if (particlesData->Get_z()[j] < 0)
                    particlesData->Get_z()[j] = particlesData->Get_z()[j] + 2 * M_PI;
            particlesData->Get_pz()[j] =  Data[1][j];

            particlesData->q[j] = this->DistribParams[0] / (this->NumbersParams[2] *
    this->DistribParams[3]);

            particlesDataZ0->Get_r()[j] = std::abs(Data[2][j]);
            particlesDataZ0->Get_z()[j] = 0;
    };*/

};
template <class PointType>
void EmitterDevice2daxs<PointType>::GenerateParticles(
    std::vector<int> indesexPerThread, int thread, int numThreads,
    std::vector<unsigned int>&                        EmptyPlaces,
    const std::shared_ptr<Particles3dcil<PointType>>& particlesData, PointType restMass,
    double charge, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData2daxs<PointType>>& grid, int flagLocate, int flagDistrib)
{
    int    nParticlesEmitter = this->NumbersParams[2];
    double energyAverage     = this->DistribParams[0];
    int    EmitPeriod        = this->NumbersParams[0];

    int perThread = nParticlesEmitter / numThreads;

    PointType nParticlesEnergyLoc;

    int nParticlesRZ     = 1;
    int nParticlesPhiR   = 1;
    int nParticlesEnergy = this->NumbersParams[1];

    if (flagDistrib == 0)
        nParticlesEnergyLoc = 1;
    else
        nParticlesEnergyLoc = nParticlesEnergy;

    if (stepNumber % EmitPeriod)
        return;

    dt = EmitPeriod * dt;
    if (flagClear == 0)
    {
        particlesData->clear();
    };

    int empty          = int(EmptyPlaces.size());
    int totalParticles = indesexPerThread.size() * nParticlesEnergyLoc;
    int nowParticles   = particlesData->NParticles();
    if (empty < totalParticles)
        particlesData->resize(totalParticles + nowParticles - empty);

    PointType dL = this->particleSource->length() / nParticlesEmitter;
    PointType L  = 0;

    //PointType EnergyAv = std::abs(energyAverage);
    //double    sigma    = sqrt(EnergyAv * std::abs(charge) / restMass);

    this->particleSource->resetSearch();

    double    vNorm;
    double    vTang;
    double    vPhi;
    double    r, z, r1, z1;
    int       cellNumb;
    double    alphaEdge;
    PointType currentFrom_dl;
    int       k    = 0;
    int       k1   = 0;
    PointType curr = 0;

    PointType    currTot = 0;
    unsigned int index;
	PointType particleData[10];
    for (int ii1 = 0; ii1 < indesexPerThread.size(); ii1++)
    {
        int    i1 = indesexPerThread[ii1];
        double L  = i1 * dL;

        this->particleSource->GetParticleOptimized(L, L + dL, 1, particleData);

        if (particleData[2] == 0)
            continue;

        alphaEdge = particleData[3];

        r = particleData[0];
        z = particleData[1];

        double alpha = particleData[4];
        cellNumb = int(particleData[6]);

        /*seachIntersictionP2.x = r1;
        seachIntersictionP2.y = z1;
        seachIntersictionP2.z = 0;

        seachIntersictionP1.x = r;
        seachIntersictionP1.y = z;
        seachIntersictionP1.z = 0;*/
	    currentFrom_dl = particleData[2];
        for (int i0 = 0; i0 < nParticlesEnergyLoc; i0++)
        {
//            vNorm = std::abs(sigma * distribution(generator));
//            vTang = sigma * distribution(generator);
//            vPhi  = sigma * distribution(generator);
			//double energy = std::abs( EmitterDeviceBase<PointType>::get_energy_distribution() ) * 1.602e-19;
			double v = 0.1;//sqrt( 2 * energy/ restMass );
            //double v    = sqrt(vNorm * vNorm + vTang * vTang + vPhi * vPhi);
            double beta = v / commtools::LIGHT_VELOCITY();
            double gamma = 1 / sqrt(1 - beta * beta);
            if (k < empty)
            {
                index          = EmptyPlaces[k];
                EmptyPlaces[k] = -1;
            }
            else
            {
                index = nowParticles + k1;
                k1++;
            }

            particlesData->Get_pr()[index] =
                gamma * (v*particleData[7]) / commtools::LIGHT_VELOCITY();
            particlesData->Get_pz()[index] =
                gamma * (v*particleData[8]) / commtools::LIGHT_VELOCITY();

            particlesData->Get_pphi()[index] = (v * particleData[9] / commtools::LIGHT_VELOCITY()) * gamma * r;
            particlesData->Get_phi()[index]  = 0;

            particlesData->Get_r()[index] = r;
            particlesData->Get_z()[index] = z;

            if (flagDistrib == 0)
            {
                particlesData->Get_pz()[index]   = 0;
                particlesData->Get_pphi()[index] = 0;
            };
            // particlesData->gamma[index] = sqrt(1 + particlesData->Get_pr()[index] *
            // particlesData->Get_pr()[index] + particlesData->Get_pz()[index] *
            // particlesData->Get_pz()[index] + (particlesData->Get_pphi()[index] /
            // particlesData->Get_r()[index])*(particlesData->Get_pphi()[index] /
            // particlesData->Get_r()[index]));

            if (currentFrom_dl < 0)
                currentFrom_dl = 0;

            currTot = currTot + currentFrom_dl / nParticlesEnergyLoc;

            particlesData->q[index] =
                currentFrom_dl * dt * Dmath::sign(charge) / nParticlesEnergyLoc;

            particlesData->cellsNumbers[index] = 0;
            particlesData->flagEmitted[index]  = 0;

            k++;
        }
    }

    if (empty < totalParticles)
        particlesData->resize(k1 + nowParticles);

    std::vector<unsigned int> remove;
    for (int i = k; i < EmptyPlaces.size(); i++)
    {
        if (EmptyPlaces[i] != -1)
        {
            remove.push_back(EmptyPlaces[i]);
            particlesData->cellsNumbers[EmptyPlaces[i]] = -1;
        }
    }

    EmptyPlaces = remove;

    //	if (std::abs(cur) / totalParticles<particlesData->avCharge)
    particlesData->avCharge = std::abs(currTot * dt) / totalParticles;
};

template <class PointType>
double EmitterDevice2daxs<PointType>::GetEmissionCurrent()
{
    return this->particleSource->GetEmissionCurrent(1);
};
template <class PointType>
std::vector<std::vector<float>> EmitterDevice2daxs<PointType>::GetCurrentDensityDistribution()
{
    return this->particleSource->GetCurrentDensityDistribution();
};
template <class PointType>
EmitterDevice2daxs<PointType>::EmitterDevice2daxs(int DistributionStyleIn)
    : EmitterDeviceBase<PointType>(DistributionStyleIn){

      };
template <class PointType>
void EmitterDevice2daxs<PointType>::GenerateParticlesLinac(
    int flagTest, int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
    const std::shared_ptr<Particles3dcil<PointType>>& particlesData, PointType restMass,
    short chargeSign, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData2daxs<PointType>>& grid, int flagLocate, int flagDistrib)
{
    /*	double energyAverage = this->DistribParams[1];
            double frequency = this->DistribParams[3];
            int particlesPerBunch = this->NumbersParams[2];

            double restEnergy = -restMass*commtools::LIGHT_VELOCITY()*commtools::LIGHT_VELOCITY() /
       commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

            double gamma = (restEnergy + std::abs(energyAverage)) / restEnergy;
            double beta = sqrt(gamma * gamma - 1) / gamma;

            if (dt == 0 || dt == 1)
            {
                    double restEnergy =
       -restMass*commtools::LIGHT_VELOCITY()*commtools::LIGHT_VELOCITY() /
       commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

                    particlesData->Get_r().resize(1);
                    particlesData->Get_z()[0] = 0;
                    particlesData->Get_pz()[0] = sqrt(pow(1 + energyAverage / (restEnergy), 2) - 1);
                    particlesData->Get_r()[0] = 0.000001;
                    particlesData->Get_pr()[0] = 0;
                    particlesData->Get_phi()[0] = 0;
                    particlesData->Get_pphi()[0] = 0;

                    return;
            }
            double t1 = particlesData->Time / commtools::LIGHT_VELOCITY();
            double t2 = particlesData->Time / commtools::LIGHT_VELOCITY() + dt;

            double z1 = 2 * M_PI*t1*frequency;
            double z2 = 2 * M_PI* t2*frequency;

            int nBunches = this->NumbersParams[3];

            if (z2 < nBunches* 2 * M_PI )
            {
                    while (z1 > 2 * M_PI)
                            z1 = z1 - 2 * M_PI;

                    while (z2 > 2 * M_PI)
                            z2 = z2 - 2 * M_PI;
            }
            //double z1 = fmod(t1*frequency, 2 * M_PI) - M_PI;
            //double z2 = fmod(t2*frequency, 2 * M_PI) - M_PI;
            std::vector<int>sliceIndexes;
            int i = 0;
            while (Data[4][i]>z2)
            {
                    i++;
                    if (i >= Data[4].size()-1)
                            break;
            }
            while (Data[4]()[i]>z1 && Data[4][i] < z2)
            {
                    sliceIndexes.push_back(i);
                    i++;
                    if (i >= Data[4]().size() - 1)
                            break;
            }

            int s = sliceIndexes.size() / nBunches;
            int rem = sliceIndexes.size() - s;

            for (int i = 0; i < rem; i++)
            {
                    int r = rand() % sliceIndexes.size();
                    sliceIndexes.erase(sliceIndexes.begin() + r);
            };

            PointType phMin = 100000;
            PointType phMax = -100000;

            for (int i = 0; i < sliceIndexes.size(); i++)
            {
                    if (Data[4][sliceIndexes[i]]>phMax)
                            phMax = Data[4][sliceIndexes[i]];

                    if (Data[4][sliceIndexes[i]]<phMin)
                            phMin = Data[4][sliceIndexes[i]];
            };





            //if (Data[4][particlesPerBunch]>z1)//there is no bunch coming at this time, just remove
       Emptyplaces
            //{//emitpiriod?
            //}

            int perThread = sliceIndexes.size() / numThreads;

            std::vector<int> sliceIndexesParallel;

            int j1 = thread*perThread;
            int j2 = (thread+1)*perThread;
            if (thread == numThreads-1)
            {
                    j2 = sliceIndexes.size();
            };

            sliceIndexesParallel.resize(j2-j1);

            for (int j = j1; j < j2; j++)
            {
                    sliceIndexesParallel[j - j1] = sliceIndexes[j];
            };


            int empty = int(EmptyPlaces.size());
            int totalParticles = sliceIndexesParallel.size();
            int nowParticles = particlesData->NParticles();
            if (empty < totalParticles)
                            particlesData->resize(totalParticles + nowParticles - empty);

            int index;

            int k1 = 0;

            for (int i = 0; i < totalParticles; i++)
            {

                    if (i < empty)
                    {
                            index = EmptyPlaces[i];
                            EmptyPlaces[i] = -1;
                    }
                    else
                    {
                            index = nowParticles + k1;
                            k1++;
                    }
                    if (phMax == phMin)
                            particlesData->Get_z()[index] = 0;
                    else
                            particlesData->Get_z()[index] =(phMax-Data[4][sliceIndexesParallel[i]])
       / (phMax - phMin)*beta*commtools::LIGHT_VELOCITY()*dt;







                    particlesData->Get_pz()[index] =
       newParticles->Get_pz()[sliceIndexesParallel[i]]; particlesData->Get_r()[index] =
       newParticles->Get_r()[sliceIndexesParallel[i]]; particlesData->Get_pr()[index] =
       newParticles->Get_pr()[sliceIndexesParallel[i]]; particlesData->Get_phi()[index] =
       newParticles->Get_phi()[sliceIndexesParallel[i]]; particlesData->Get_pphi()[index] =
       newParticles->Get_pphi()[sliceIndexesParallel[i]]; particlesData->q[index] =
       newParticles->q[sliceIndexesParallel[i]];
            //	particlesData->gamma[index] = newParticles->gamma[sliceIndexesParallel[i]];
            }

            std::vector<unsigned int> remove;
            for (int i = totalParticles; i < EmptyPlaces.size(); i++)
            {
                    if (EmptyPlaces[i] != -1)
                    {
                            remove.push_back(EmptyPlaces[i]);
                            particlesData->cellsNumbers[EmptyPlaces[i]] = -1;
                    }
            }


            EmptyPlaces = remove;*/
}

template <class PointType>
void EmitterDevice2daxs<PointType>::GenerateSyncParticle(
    const std::shared_ptr<Particles3dcil<PointType>>& particlesData, PointType restMass)
{
    double restEnergy = -restMass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    particlesData->resize(1);
    particlesData->Get_z()[0]    = 0;
    particlesData->Get_pz()[0]   = sqrt(pow(1 + this->DistribParams[1] / (restEnergy), 2) - 1);
    particlesData->Get_r()[0]    = 0.000001;
    particlesData->Get_pr()[0]   = 0;
    particlesData->Get_phi()[0]  = 0;
    particlesData->Get_pphi()[0] = 0;
};

template <class PointType>
std::vector<double> EmitterDevice2daxs<PointType>::GetAdditionalSourceInf()
{
    std::vector<double> res;
    return res;
}
template <class PointType>
void EmitterDevice2daxs<PointType>::SetAdditionalSourceInf(std::vector<double> inf){};
template <class PointType>
void EmitterDevice2daxs<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData2daxs<PointType>>&             grid)
{
    this->boundaryList     = in;
    this->emitterInitParam = parametersIn;
    this->particleSource->InitEmissionBoundary(boundaryIn, grid, parametersIn, error);
};
