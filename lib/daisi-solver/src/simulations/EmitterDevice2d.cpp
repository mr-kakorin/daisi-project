#include "EmitterDevice2d.h"
#include "Dmath.h"
#include "EmitterDevice2daxsVectorized.h"
#include "Geom.h"
#include "GridData.h"
#include "Particle.h"
#include "ParticleSource.h"
#include "Tools.h"
#include <common_tools/constants.h>
#include <random>

static std::default_random_engine       generator;
static std::normal_distribution<double> distribution(0, 1);

template class EmitterDevice2d<float>;
template class EmitterDevice2d<double>;

std::vector<float> steps(512);
std::vector<float> tmp(512);
std::vector<float> tmp1(512);

std::vector<std::vector<float>> stepsT(10);
std::vector<std::vector<float>> tmpT(10);
std::vector<std::vector<float>> tmp1T(10);

double IntegrateVelocity(double a, double b, double t0, double m, double charge, int thread)
{
    int    s    = 300;
    double step = (b - a) / s;

    double x = a;
    stepsT[thread].resize(512);
    tmpT[thread].resize(512);
    stepsT[thread][0] = a;

    for (int i = 1; i < s + 1; i++)
    {
        stepsT[thread][i] = stepsT[thread][i - 1] + step;
    };

    for (int i = 0; i < s + 1; i++)
    {

        tmpT[thread][i] =
            sqrt(m / (2 * commtools::PI() * t0 * std::abs(charge))) *
            exp(-(m * stepsT[thread][i] * stepsT[thread][i]) / (2 * t0 * std::abs(charge)));
    }

    double result = 0;

    for (int i = 0; i < s; i++)
    {
        result = result + tmpT[thread][i] + tmpT[thread][i + 1];
    };
    result = result * step / 2;

    return result;
};

template <class PointType>
void EmitterDevice2d<PointType>::GenerateParticles(
    std::vector<int> indesexPerThread, int thread, int numThreads,
    std::vector<unsigned int>&                     EmptyPlaces,
    const std::shared_ptr<Particles2d<PointType>>& particlesData, PointType restMass, double charge,
    int flagClear, double dt, int stepNumber, const std::shared_ptr<GridData2d<PointType>>& grid,
    int flagLocate, int flagDistrib)
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

    PointType EnergyAv = std::abs(energyAverage);
    double    sigma    = sqrt(EnergyAv * std::abs(charge) / restMass);

    this->particleSource->resetSearch();

    double    vNorm;
    double    vTang;
    double    vPhi;
    double    x, y, x1, y1;
    int       cellNumb;
    double    alphaEdge;
    PointType currentFrom_dl;
    int       k    = 0;
    int       k1   = 0;
    PointType curr = 0;

    PointType    currTot = 0;
    unsigned int index;

    for (int ii1 = 0; ii1 < indesexPerThread.size(); ii1++)
    {
        int    i1 = indesexPerThread[ii1];
        double L  = i1 * dL;

        std::vector<PointType> tmp = this->particleSource->GetParticle(L, L + dL, 0);

        if (tmp[2] == 0)
            continue;

        alphaEdge = tmp[3];

        x = tmp[0];
        y = tmp[1];

        x1 = tmp[4];
        y1 = tmp[5];

        double alpha = tmp[4];
        // double alpha = commtools::PI() - atan((y1 - y) / (x1 - x));

        cellNumb = int(tmp[6]);

        /*seachIntersictionP2.x = r1;
        seachIntersictionP2.y = z1;
        seachIntersictionP2.z = 0;

        seachIntersictionP1.x = r;
        seachIntersictionP1.y = z;
        seachIntersictionP1.z = 0;*/
        for (int i0 = 0; i0 < nParticlesEnergyLoc; i0++)
        {
            vTang = sigma * distribution(generator);
            vNorm = std::abs(sigma * distribution(generator));

            // vTang = 0;

            double v    = sqrt(vNorm * vNorm + vTang * vTang);
            double beta = v / commtools::LIGHT_VELOCITY();

            double gamma = 1 / sqrt(1 - beta * beta);

            currentFrom_dl = tmp[2];

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

            particlesData->Get_px()[index] =
                gamma * (vNorm * std::cos(alpha) - vTang * std::sin(alpha)) / commtools::LIGHT_VELOCITY();
            particlesData->Get_py()[index] =
                gamma * (vNorm * std::sin(alpha) + vTang * std::cos(alpha)) / commtools::LIGHT_VELOCITY();

            particlesData->Get_x()[index] = x;
            particlesData->Get_y()[index] = y;

            particlesData->Get_z()[index]  = 0;
            particlesData->Get_pz()[index] = 0;

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
void EmitterDevice2d<PointType>::GenerateParticles(
    std::vector<int> indesexPerThread, int thread, int numThreads,
    std::vector<unsigned int>&                          EmptyPlaces,
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesData, PointType restMass,
    double charge, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData2dpolar<PointType>>& grid, int flagLocate, int flagDistrib)
{
    int    nParticlesEmitter = this->NumbersParams[3];
    double energyAverage     = this->DistribParams[0];
    int    EmitPeriod        = this->NumbersParams[1];

    int perThread = nParticlesEmitter / numThreads;

    PointType nParticlesRZLoc;
    PointType nParticlesPhiRLoc;

    int nParticlesRZ     = 1;
    int nParticlesPhiR   = 1;
    int nParticlesEnergy = this->NumbersParams[2];

    int i1_0 = perThread * thread;
    int i1_1 = perThread * thread + perThread;

    if (thread == numThreads - 1)
        i1_1 = nParticlesEmitter;

    if (stepNumber % EmitPeriod)
    {
        particlesData->removeParticle(EmptyPlaces);
        return;
    }

    int nParticlesEnergyLoc;
    int nParticlesXYLoc;

    if (flagDistrib == 0)
    {
        nParticlesEnergyLoc = 1;
        nParticlesXYLoc     = 1;
    }
    else
    {
        nParticlesEnergyLoc = nParticlesEnergy;
        nParticlesXYLoc     = nParticlesXY;
    };

    dt = EmitPeriod * dt;
    if (flagClear == 0)
    {
        particlesData->clear();
    };
    int empty = int(EmptyPlaces.size());

    int totalParticles = indesexPerThread.size() * nParticlesEnergyLoc;
    int nowParticles   = particlesData->NParticles();
    if (empty < totalParticles)
        particlesData->resize(totalParticles + nowParticles - empty);

    PointType dL = this->particleSource->length() / nParticlesEmitter;
    PointType L  = 0;

    PointType particleEnergy;
    PointType restEnergy;
    PointType gamma;
    PointType beta;
    PointType pTotal;

    PointType x;
    PointType y;
    PointType x1;
    PointType y1;
    PointType rtmp;
    PointType phitmp;
    PointType currentFrom_dl;
    PointType CphXY;
    PointType phXY0;
    PointType phXY1;
    PointType phXY;
    PointType CphRPhi;
    PointType phPhiR;
    PointType current;
    PointType dphiXY = commtools::PI() / nParticlesXYLoc;

    int       k    = 0;
    int       k1   = 0;
    PointType curr = 0;

    PointType    currTot = 0;
    unsigned int index;

    this->particleSource->resetSearch();

    this->particleSource->resetSearch();

    PointType              cur = 0;
    PointType              alphaEdge;
    DGeo::Point<PointType> seachIntersictionP1;
    DGeo::Point<PointType> seachIntersictionP2;
    DGeo::Edge<PointType>  EdgeIntersection;
    DGeo::Point<PointType> startPoint;
    int                    cellNumb;

    PointType EnergyAv = std::abs(energyAverage);

    PointType velocityRight = 0;
    PointType dAv           = sqrt(restMass / (2 * commtools::PI() * EnergyAv * std::abs(charge))) *
                    exp(-(restMass * 0 * 0) / (2 * EnergyAv * std::abs(charge)));

    while (1)
    {
        velocityRight = velocityRight + 1000;
        PointType dAvT =
            sqrt(restMass / (2 * commtools::PI() * EnergyAv * std::abs(charge))) *
            exp(-(restMass * velocityRight * velocityRight) / (2 * EnergyAv * std::abs(charge)));
        if (dAvT < 0.01 * dAv)
            break;
    };

    PointType velocityLeft = -velocityRight;

    PointType dvPhi = (velocityRight - velocityLeft) / nParticlesEnergyLoc;

    restEnergy = -restMass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                 commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    PointType gammaR  = (restEnergy + std::abs(energyAverage)) / restEnergy;
    PointType betaR   = sqrt(gammaR * gammaR - 1) / gammaR;
    PointType pTotalR = betaR * gammaR;

    PointType vPhi;
    PointType vR;
    PointType v;

    gammaR  = (restEnergy + std::abs(energyAverage)) / restEnergy;
    betaR   = sqrt(gammaR * gammaR - 1) / gammaR;
    pTotalR = betaR * gammaR;

    for (int ii1 = 0; ii1 < indesexPerThread.size(); ii1++)
    {
        int i1 = indesexPerThread[ii1];
        L      = i1 * dL;

        std::vector<PointType> tmp = this->particleSource->GetParticle(L, L + dL, 0);

        for (int i0 = 0; i0 < nParticlesEnergyLoc; i0++)
        {

            int sign = energyAverage / std::abs(energyAverage);

            PointType c1 =
                IntegrateVelocity(velocityLeft + i0 * dvPhi, velocityLeft + (i0 + 1) * dvPhi,
                                  EnergyAv, restMass, charge, thread);

            if (nParticlesEnergyLoc == 1)
            {
                c1   = 1;
                vPhi = 0;
            }
            else
            {
                vPhi = velocityLeft + i0 * dvPhi + dvPhi / 2;
            }

            gamma = (restEnergy + std::abs(EnergyAv)) / restEnergy;
            beta  = sqrt(gamma * gamma - 1) / gamma;
            vR    = beta * commtools::LIGHT_VELOCITY();

            v    = sqrt(vR * vR + vPhi * vPhi);
            beta = v / commtools::LIGHT_VELOCITY();

            gamma = 1 / sqrt(1 - beta * beta);

            pTotal = beta * gamma;

            x              = tmp[0];
            y              = tmp[1];
            currentFrom_dl = tmp[2];

            alphaEdge = tmp[3];

            x = tmp[0];
            y = tmp[1];

            x1       = tmp[4];
            y1       = tmp[5];
            cellNumb = int(tmp[6]);

            seachIntersictionP2.x = x1;
            seachIntersictionP2.y = y1;
            seachIntersictionP2.z = 0;

            seachIntersictionP1.x = x;
            seachIntersictionP1.y = y;
            seachIntersictionP1.z = 0;

            EdgeIntersection.point1 = seachIntersictionP1;
            EdgeIntersection.point2 = seachIntersictionP2;

            PointType alphaN = EdgeIntersection.alpha();

            curr = curr + currentFrom_dl;
            // CphXY = 0.5 * currentFrom_dl / IntegrateCurrent(0, commtools::PI() / 2, phiXY);

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

            phXY0 = -commtools::PI() / 2;

            if (flagLocate)
            {
                grid->SearchIntersectionWithEdge(cellNumb, EdgeIntersection, &startPoint);
                particlesData->GetPointerToPosition1()[index] = startPoint.x;
                particlesData->GetPointerToPosition2()[index] = startPoint.y;
            }
            else
            {
                particlesData->GetPointerToPosition1()[index] = x;
                particlesData->GetPointerToPosition2()[index] = y;
            }

            Dmath::Cartesian2Polar(particlesData->GetPointerToPosition1()[index],
                                   particlesData->GetPointerToPosition2()[index], rtmp, phitmp);

            particlesData->Get_phi()[index]               = phitmp;
            particlesData->GetPointerToPosition1()[index] = rtmp;

            particlesData->GetPointerToMomentum1()[index] =
                (vR / commtools::LIGHT_VELOCITY()) * gamma * Dmath::sign(alphaEdge);
            particlesData->Get_pphi()[index] = (vPhi / commtools::LIGHT_VELOCITY()) * gamma *
                                               particlesData->GetPointerToPosition1()[index];

            if (particlesData->GetPointerToPosition1()[index] > 0.05 &&
                particlesData->GetPointerToPosition1()[index] < 0.12)
            {
                //	if (alphaEdge < 0)
                //		alphaEdge = std::abs(2 * commtools::PI() - phitmp) - commtools::PI(); //
                //�������!! 	else 		alphaEdge = std::abs(2 * commtools::PI() - phitmp); //
                //�������!!

                phXY = commtools::PI();

                particlesData->SetCartesianMomentum(
                    index, {pTotal * std::cos(phXY + alphaEdge), pTotal * std::sin(phXY + alphaEdge)});

                //	particlesData->SetCartesianMomentum(index, -pTotalR, pTotal);
            }

            double beta = particlesData->GetBeta(index);

            // particlesData->gamma[index] = sqrt(1 + particlesData->GetPointerToMomentum1()[index]
            // * particlesData->GetPointerToMomentum1()[index] + (particlesData->Get_pphi()[index] /
            // particlesData->GetPointerToPosition1()[index])*(particlesData->Get_pphi()[index] /
            // particlesData->GetPointerToPosition1()[index]));

            currTot = currTot + currentFrom_dl * c1;

            //	if (current < 0)
            //		current = 0;
            particlesData->q[index]            = currentFrom_dl * c1 * dt * Dmath::sign(charge);
            particlesData->cellsNumbers[index] = 0;
            particlesData->flagEmitted[index]  = 0;

            //	phXY0 = phXY1;

            k++;
        };
    };

    //	particlesData->resize(k1);
    //// 0 r 1 phi 2 x 3 y 4 phiReal 5 xReal 6 yReal,

    particlesData->positions[4] = particlesData->positions[1];
    particlesData->positions[5] = particlesData->positions[2];
    particlesData->positions[6] = particlesData->positions[3];

    std::vector<unsigned int> remove;
    for (int i = k; i < EmptyPlaces.size(); i++)
        remove.push_back(EmptyPlaces[i]);

    EmptyPlaces = remove;
    // particlesData->removeParticle(remove);

    //	if (std::abs(cur) / totalParticles<particlesData->avCharge)
    particlesData->avCharge = std::abs(currTot * dt) / totalParticles;

    this->particleSource->sourceCurrent = currTot;
};

template <class PointType>
void EmitterDevice2d<PointType>::GenerateParticlesLinac(
    int flagTest, int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
    const std::shared_ptr<Particles2d<PointType>>& particlesData, PointType restMass,
    short chargeSign, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData2d<PointType>>& grid, int flagLocate, int flagDistrib)
{
    double energyAverage     = this->DistribParams[1];
    double frequency         = this->DistribParams[3];
    int    particlesPerBunch = this->NumbersParams[0];

    double restEnergy = -restMass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(energyAverage)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    std::vector<int> sliceIndexesParallel;
    double           t1 = particlesData->Time / commtools::LIGHT_VELOCITY();
    double           t2 = particlesData->Time / commtools::LIGHT_VELOCITY() + dt;

    double phMin;
    double phMax;

    this->GetSliceIndexes(sliceIndexesParallel, flagTest, t1, t2, phMin, phMax, numThreads, thread);

    int empty = int(EmptyPlaces.size());

    int totalParticles = sliceIndexesParallel.size();

    int nowParticles = particlesData->NParticles();
    if (empty < totalParticles)
        particlesData->resize(totalParticles + nowParticles - empty);

    int k1    = 0;
    int k     = 0;
    int index = -1;
    for (k = 0; k < totalParticles; k++)
    {

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

        if (phMax == phMin)
            particlesData->Get_z()[index] = 0;
        else
        {
            double frequency = this->DistribParams[3];
            double dtBunch   = (phMax - phMin) / (2 * commtools::PI() * frequency);
            double add       = beta * commtools::LIGHT_VELOCITY() *
                         (2 * commtools::PI() * t2 * frequency - phMax) /
                         (2 * commtools::PI() * frequency);
            particlesData->Get_z()[index] =
                add +
                ((phMax - (this->Data)[4][sliceIndexesParallel[k]]) / (phMax - phMin)) * beta *
                    commtools::LIGHT_VELOCITY() * dtBunch;
        }

        particlesData->Get_x()[index] = (this->Data)[0][sliceIndexesParallel[k]];
        particlesData->Get_y()[index] = (this->Data)[2][sliceIndexesParallel[k]];

        particlesData->Get_px()[index] = (this->Data)[1][sliceIndexesParallel[k]] * beta;
        particlesData->Get_py()[index] = (this->Data)[3][sliceIndexesParallel[k]] * beta;
        particlesData->Get_pz()[index] = (this->Data)[5][sliceIndexesParallel[k]];

        particlesData->q[index] =
            this->DistribParams[0] / (this->NumbersParams[0] * commtools::LIGHT_VELOCITY() * beta);
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
    if (-1 != index)
        particlesData->avCharge = particlesData->q[index];
}
template <class PointType>
void EmitterDevice2d<PointType>::GenerateParticlesLinac(
    int flagTest, int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesData, PointType restMass,
    short chargeSign, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData2dpolar<PointType>>& grid, int flagLocate, int flagDistrib)
{
    int temp1 = 2;
}
template <class PointType>
double EmitterDevice2d<PointType>::GetEmissionCurrent()
{
    return this->particleSource->GetEmissionCurrent(0);
};

template <class PointType>
EmitterDevice2d<PointType>::EmitterDevice2d(int DistributionStyleIn)
    : EmitterDeviceBase<PointType>(DistributionStyleIn){};

template <class PointType>
EmitterDevice2d<PointType>::EmitterDevice2d(){};
template <class PointType>
std::vector<std::vector<float>> EmitterDevice2d<PointType>::GetCurrentDensityDistribution()
{
    return this->particleSource->GetCurrentDensityDistribution();
};

template <class PointType>
void EmitterDevice2d<PointType>::GenerateSyncParticle(
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesData, PointType restMass){

};
template <class PointType>
void EmitterDevice2d<PointType>::GenerateSyncParticle(
    const std::shared_ptr<Particles2d<PointType>>& particlesData, PointType restMass){

};
template <class PointType>
std::vector<double> EmitterDevice2d<PointType>::GetAdditionalSourceInf()
{
    std::vector<double> res;
    return res;
}
template <class PointType>
void EmitterDevice2d<PointType>::SetAdditionalSourceInf(std::vector<double> inf){};

template <class PointType>
void EmitterDevice2d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaryIn,
    const std::shared_ptr<GridData2d<PointType>>&                 grid)
{
    this->boundaryList     = in;
    this->emitterInitParam = parametersIn;
    this->particleSource->InitEmissionBoundary(boundaryIn, grid, parametersIn, error);
};
template <class PointType>
void EmitterDevice2d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaryIn,
    const std::shared_ptr<GridData2dpolar<PointType>>&            grid)
{
    this->boundaryList     = in;
    this->emitterInitParam = parametersIn;
    this->particleSource->InitEmissionBoundary(boundaryIn, grid, parametersIn, error);
};
template <class PointType>
void EmitterDevice2d<PointType>::PreliminaryGeneration(
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesDataZ0, PointType restMass){

};
template <class PointType>
void EmitterDevice2d<PointType>::PreliminaryGeneration(
    const std::shared_ptr<Particles2d<PointType>>& particlesDataZ0, PointType restMass)
{
    this->GenerateEllipses(restMass);

    int    particlesPerBunch = this->NumbersParams[0] * this->NumbersParams[1];
    double energyAverage     = this->DistribParams[1];
    double restEnergy = -restMass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(energyAverage)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    particlesDataZ0->resize(particlesPerBunch);

    std::vector<double> en(particlesPerBunch);
    for (int j = 0; j < particlesPerBunch; j++)
    {
        /*particlesData->Get_x()[j] = (this->Data)[0][j];
        particlesData->Get_y()[j] = (this->Data)[2][j];
        particlesData->Get_z()[j] = (this->Data)[4][j];

        particlesData->Get_px()[j] = (this->Data)[1][j] * beta;
        particlesData->Get_py()[j] = (this->Data)[3][j] * beta;
        particlesData->Get_pz()[j] = (this->Data)[5][j];

        particlesData->q[j] = this->DistribParams[0] / (this->NumbersParams[0] *
        commtools::LIGHT_VELOCITY() * beta);*/

        particlesDataZ0->Get_x()[j] = (this->Data)[0][j];
        particlesDataZ0->Get_y()[j] = (this->Data)[2][j];
        particlesDataZ0->Get_z()[j] = 0;
    };
};

template void
EmitterDevice2d<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);

template void EmitterDevice2d<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void
EmitterDevice2d<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);

template void EmitterDevice2d<float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void EmitterDevice2d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& nParticlesXY;
    ar& phiXY;
    ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
}
template <class PointType>
template <class Archive>
void EmitterDevice2d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& nParticlesXY;
    ar& phiXY;
    ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
}