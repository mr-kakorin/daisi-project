#include "EmitterDevice3d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "EmitterDevice2daxsVectorized.h"
#include "Geom.h"
#include "GridData.h"
#include "Particle.h"
#include "ParticleSource.h"
#include "Tools.h"
#include <Constants.h>


template <class PointType>
void EmitterDevice3d<PointType>::GenerateParticles(
    std::vector<int> indesexPerThread, int thread, int numThreads,
    std::vector<unsigned int>&                     EmptyPlaces,
    const std::shared_ptr<Particles3d<PointType>>& particlesData, PointType restMass, double charge,
    int flagClear, double dt, int stepNumber, const std::shared_ptr<GridData3d<PointType>>& grid,
    int flagLocate, int flagDistrib){

    /*
    double energyAverage = DistribParams[0];
    int EmitPeriod

    int perThread = nParticlesZ / numThreads;

    int i1_0 = perThread*thread;
    int	i1_1 = perThread*thread + perThread;

    if (thread == numThreads - 1)
            i1_1 = nParticlesZ;


    if (stepNumber%EmitPeriod)
    {
            particlesData->removeParticle(EmptyPlaces);
            return;
    }

    dt = EmitPeriod*dt;
    if (flagClear == 0)
    {
            particlesData->clear();
    };

    int nParticlesEnergyLoc;
    int nParticlesXYLoc;


    if (flagDistrib == 0)
    {
            nParticlesEnergyLoc = 1;
            nParticlesXYLoc = 1;
    }
    else
    {
            nParticlesEnergyLoc = nParticlesEnergy;
            nParticlesXYLoc = nParticlesXY;
    };

    int empty = int(EmptyPlaces.size());
    int totalParticles = (i1_1 - i1_0)*nParticlesEmitter*nParticlesEnergyLoc*nParticlesXYLoc;
    int nowParticles = particlesData->NParticles();
    if (empty<totalParticles)
            particlesData->resize()(totalParticles + nowParticles - empty);

    PointType dL = this->particleSources[0]->length() / nParticlesEmitter;
    PointType L = 0;


    PointType particleEnergy;
    PointType restEnergy;
    PointType gamma;
    PointType beta;
    PointType pTotal;
    PointType x;
    PointType y;
    PointType x1;
    PointType y1;
    PointType currentFrom_dl;
    PointType CphXY;
    PointType phXY0;
    PointType phXY1;
    PointType phXY;
    PointType CphRPhi;
    PointType phPhiR;
    PointType current;
    PointType dphiXY = PI() / nParticlesXYLoc;
    PointType dZ = (Z2-Z1) / nParticlesZ;
    PointType z = Z1 + dZ/2;

    int k = 0;
    int k1 = 0;
    PointType curr = 0;

    PointType currTot = 0;
    unsigned int index;



    PointType cur = 0;
    PointType alphaEdge;
    DGeo::Point<PointType> seachIntersictionP1;
    DGeo::Point<PointType> seachIntersictionP2;
    DGeo::Edge<PointType> EdgeIntersection;
    DGeo::Point<PointType> startPoint;
    int cellNumb;


    for (int i00 = i1_0; i00 < i1_1; i00++)
    {
            int sourceNumber;
            for (sourceNumber = 0; sourceNumber < zArray.size(); sourceNumber++)
            {
                    if (z >= zArray[sourceNumber] && z <= zArray[sourceNumber + 1])
                            break;
            };

            this->particleSources[sourceNumber].resetSearch();
            this->particleSources[sourceNumber+1].resetSearch();

            for (int i0 = 0; i0 < nParticlesEnergyLoc; i0++)
            {
                    int sign = energyAverage / std::abs(energyAverage);
                    particleEnergy = std::abs(energyAverage); // ��� ����� ������ �������� ����������
    �������������, �� ���-������ � ������ ��� restEnergy =
    -restMass*LIGHT_VELOCITY()*LIGHT_VELOCITY() /
    ELECTRON_CHARGE(); // ������� ����� ������� � ��������������� gamma = (restEnergy +
    std::abs(energyAverage)) / restEnergy; beta = sqrt(gamma * gamma - 1) / gamma; pTotal = beta * gamma;
    for (int i1 = 0; i1 < nParticlesEmitter; i1++)
                    {
                            L = i1 * dL;

                            std::vector<PointType> tmp =
    this->particleSources[sourceNumber].GetParticle(L, L + dL, 0); std::vector<PointType> tmp1 =
    this->particleSources[sourceNumber].GetParticle(L, L + dL, 0);

                            x = tmp[0];
                            y = tmp[1];

                            double w1 = (zArray[sourceNumber + 1] - z) / (zArray[sourceNumber + 1] -
    zArray[sourceNumber]); double w2 = 1 - w1;

                            currentFrom_dl = tmp[2] * w1 + tmp1[2]*w2;
                            alphaEdge = tmp[3];

                            alphaEdge = tmp[3];

                            x = tmp[0];
                            y = tmp[1];

                            x1 = tmp[4];
                            y1 = tmp[5];
                            cellNumb = int(tmp[6]);

                            seachIntersictionP2.x = x1;
                            seachIntersictionP2.y = y1;
                            seachIntersictionP2.z = 0;

                            seachIntersictionP1.x = x;
                            seachIntersictionP1.y = y;
                            seachIntersictionP1.z = 0;


                            curr = curr + currentFrom_dl;
                            CphXY = 0.5 * currentFrom_dl / IntegrateCurrent(0, PI() / 2,
    phiXY);

                            phXY0 = -PI() / 2;
                            for (int i2 = 0; i2 < nParticlesXYLoc; i2++)
                            {
                                    phXY1 = phXY0 + dphiXY;
                                    phXY = (phXY1 + phXY0) / 2;

                                    current = CphXY * IntegrateCurrent(phXY0, phXY1, phiXY);

                                    if (k < empty)
                                            index = EmptyPlaces[k];
                                    else
                                    {
                                            index = nowParticles + k1;
                                            k1++;
                                    }



                                    particlesData->GetPointerToMomentum1()[index] = pTotal*std::cos(phXY
    + alphaEdge); particlesData->GetPointerToMomentum2()[index] = pTotal*std::sin(phXY + alphaEdge);
                                    particlesData->GetPointerToMomentum3()[index] = 0;

                                    seachIntersictionP2 = seachIntersictionP2.rotate(phXY,
    seachIntersictionP1);

                                    EdgeIntersection.point1 = seachIntersictionP1;
                                    EdgeIntersection.point2 = seachIntersictionP2;



                                    if (flagLocate)
                                    {
                                            grid->SearchIntersectionWithEdge(cellNumb,
    EdgeIntersection, &startPoint); particlesData->GetPointerToPosition1()[index] = startPoint.x;
                                            particlesData->GetPointerToPosition2()[index] =
    startPoint.y;
                                    }
                                    else
                                    {
                                            particlesData->GetPointerToPosition1()[index] = x;
                                            particlesData->GetPointerToPosition2()[index] = y;
                                    }

                                    particlesData->GetPointerToPosition2()[index] = z;

                                    currTot = currTot + current;
                                    if (current < 0)
                                            current = 0;
                                    particlesData->q[index] = current*dt*Dmath::sign(charge);
                                    particlesData->cellsNumbers[index] = 0;
                                    particlesData->flagEmitted[index] = 0;

                                    particlesData->gamma[index] = sqrt(1 +
    particlesData->GetPointerToMomentum1()[index]
    * particlesData->GetPointerToMomentum1()[index] + particlesData->GetPointerToMomentum2()[index]
    * particlesData->GetPointerToMomentum2()[index] + particlesData->GetPointerToMomentum3()[index]
    * particlesData->GetPointerToMomentum3()[index]);


                                    phXY0 = phXY1;

                                    k++;


                            };
                    };
            };
            z = z + dZ;

    }
    std::vector<unsigned int> remove;
    for (int i = k; i < EmptyPlaces.size(); i++)
            remove.push_back(EmptyPlaces[i]);


    EmptyPlaces = remove;
    //particlesData->removeParticle(remove);


    //	if (std::abs(cur) / totalParticles<particlesData->avCharge)
    particlesData->avCharge = std::abs(currTot*dt) / totalParticles;*/
}

template <class PointType>
void EmitterDevice3d<PointType>::GenerateParticlesLinac(
    int flagTest, int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
    const std::shared_ptr<Particles3d<PointType>>& particlesData, PointType restMass,
    short chargeSign, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridData3d<PointType>>& grid, int flagLocate, int flagDistrib){
    /*	double energyAverage = DistribParams[0];

            int totalParticles = NumbersParams[0] / numThreads; // ������� ������ ��������� ��
       ������ ���� if (thread == numThreads - 1){ totalParticles = NumbersParams[0] -
       totalParticles*(numThreads - 1);
            }

            DGeo::Point<PointType> begin, end, delta;
            double radius = 0.005;

            double totalCurrent = DistribParams[1];

            double restEnergy = -restMass*LIGHT_VELOCITY()*LIGHT_VELOCITY() /
       ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

            int nowParticles = particlesData->NParticles(); // ������� �� ������
            int empty = int(EmptyPlaces.size());

            if (empty<totalParticles)
                    particlesData->resize(totalParticles + nowParticles - empty); // �����������
       ������ ��� ����� �������, ���� ��� �����

            double currTot = totalCurrent;
            int k1 = 0;
            int k = 0;
            int index;

            for (int i = 0; i < totalParticles; i++)
            {
                    int flag = 0;
                    while (flag == 0){
                            begin.x = (rand() % 1000) / 10000.0;
                            begin.y = (rand() % 1000) / 10000.0;
                            if (sqrt(begin.y*begin.y + begin.x*begin.x) < radius){
                                    flag = 1;
                            }
                    }
                    begin.z = 0.0;

                    end.x = begin.x;
                    end.y = begin.y;
                    end.z = 0.66;

                    delta.x = end.x - begin.x;
                    delta.y = end.y - begin.y;
                    delta.z = end.z - begin.z;


                    double q = totalCurrent*dt; // ����� �����, ������� ������ ������ - ��� ��� ��
       �����


                    if (k < empty)
                    {
                            index = EmptyPlaces[k];
                            EmptyPlaces[k] = -1;
                    }
                    else
                    {
                            index = nowParticles + k1;
                            k1++;
                    }

                    //
                    particlesData->GetPointerToPosition1()[index] = begin.x;
                    particlesData->GetPointerToPosition2()[index] = begin.y;
                    particlesData->GetPointerToPosition2()[index] = begin.z;
                    particlesData->q[index] = q / totalParticles; // � ���� �� ��� �����������
       ���������� �����

                    // ��������, ��� � ���� ������ ���������� ������� , � ��� ��� � pz.
                    //��� ������, ��������, ����� ������� ��������� � ������� ��������� ���������
       ��������� �������� vx, vy, vz. double gamm� = (restEnergy + std::abs(energyAverage)) /
       restEnergy;// ������ ������� double beta = sqrt(gamm�
       * gamm� - 1) / gamm�;// ����������� �������� double pTotal = beta * gamm�;//�����������
       �������

                    // ���� �������� ����� �������
                    double n_ = pTotal / sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);

                    particlesData->GetPointerToMomentum1()[index] = delta.x*n_;
                    particlesData->GetPointerToMomentum2()[index] = delta.y*n_;
                    particlesData->GetPointerToMomentum3()[index] = delta.z*n_;

                    // ��� ����� ������� �����������

                    //particlesData->gamma[index] = gamm�;

                    particlesData->cellsNumbers[index] = 0;
                    particlesData->flagEmitted[index] = 0;
                    k++;
            }

            //��������, �������, ������� ������ ����� ���� ������� �����
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

            particlesData->avCharge = std::abs(currTot*dt) / totalParticles;*/

}

template <class PointType>
void EmitterDevice3d<PointType>::PreliminaryGeneration(
    const std::shared_ptr<Particles3d<PointType>>& particlesDataZ0, PointType restMass){
    /*double energyAverage = DistribParams[0];

    int totalParticles =10 * NumbersParams[0];


    DGeo::Point<PointType> begin, end, delta;
    double radius = 0.005;

    double totalCurrent = DistribParams[1];

    double restEnergy = -restMass*LIGHT_VELOCITY()*LIGHT_VELOCITY() /
    ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    int nowParticles = particlesData->NParticles(); // ������� �� ������

    double currTot = totalCurrent;
    int k1 = 0;
    int k = 0;
    int index;

    for (int i = 0; i < totalParticles; i++)
    {
            int flag = 0;
            while (flag == 0){
                    begin.x = (rand() % 1000) / 10000.0;
                    begin.y = (rand() % 1000) / 10000.0;
                    if (sqrt(begin.y*begin.y + begin.x*begin.x) < radius){
                            flag = 1;
                    }
            }
            begin.z = 0.0;

            end.x = begin.x;
            end.y = begin.y;
            end.z = 0.66;

            delta.x = end.x - begin.x;
            delta.y = end.y - begin.y;
            delta.z = end.z - begin.z;



            index = i;


            //
            particlesData->GetPointerToPosition1()[index] = begin.x;
            particlesData->GetPointerToPosition2()[index] = begin.y;
            particlesData->GetPointerToPosition2()[index] = begin.z;

            // ��������, ��� � ���� ������ ���������� ������� , � ��� ��� � pz.
            //��� ������, ��������, ����� ������� ��������� � ������� ��������� ��������� ���������
    �������� vx, vy, vz. double gamm� = (restEnergy + std::abs(energyAverage)) / restEnergy;// ������
    ������� double beta = sqrt(gamm� * gamm� - 1) / gamm�;// ����������� �������� double pTotal =
    beta * gamm�;//����������� �������

            // ���� �������� ����� �������
            double n_ = pTotal / sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);

            particlesData->GetPointerToMomentum1()[index] = delta.x*n_;
            particlesData->GetPointerToMomentum2()[index] = delta.y*n_;
            particlesData->GetPointerToMomentum3()[index] = delta.z*n_;

            // ��� ����� ������� �����������

    //	particlesData->gamma[index] = gamm�;

            particlesData->cellsNumbers[index] = 0;
            particlesData->flagEmitted[index] = 0;
            k++;
    }*/

}

template <class PointType>
EmitterDevice3d<PointType>::EmitterDevice3d(int DistributionStyleIn)
    : EmitterDeviceBase<PointType>(DistributionStyleIn)
{
    nParticlesXY = 1;
    phiXY        = 0.0;
    Z1           = 0;
    Z2           = 0;
    nParticlesZ  = 1;
    ndZ          = 1;
}

template <class PointType>
std::vector<std::vector<float>> EmitterDevice3d<PointType>::GetCurrentDensityDistribution()
{
    return this->particleSource->GetCurrentDensityDistribution();
}

template <class PointType>
double EmitterDevice3d<PointType>::getErAverage()
{
    double res = 0;

    if (this->particleSources.size() == 0)
        return 0;

    for (int i = 0; i < this->particleSources.size(); i++)
        res = res + this->particleSources[i]->getErAverage();

    return res / double(this->particleSources.size());
}

template <class PointType>
double EmitterDevice3d<PointType>::GetEmissionCurrent()
{
    double res = 0;

    if (this->particleSources.size() == 0)
        return 0;

    for (int i = 0; i < this->particleSources.size(); i++)
        res = res + this->particleSources[i]->GetEmissionCurrent(0);

    return res;
}

template <class PointType>
double EmitterDevice3d<PointType>::GetSourceSize()
{
    if (this->particleSources.size() == 0)
        return 0;
    return this->particleSources[0]->sourceSurface.back().curveLength;
}

template <class PointType>
std::vector<std::shared_ptr<ParticleSource2d<PointType>>>&
EmitterDevice3d<PointType>::GetParticleSource()
{
    return this->particleSources;
}

template <class PointType>
std::vector<std::shared_ptr<ParticleSource2d<PointType>>>&
EmitterDevice3d<PointType>::GetParticleSources()
{
    return this->particleSources;
}

template <class PointType>
std::vector<double> EmitterDevice3d<PointType>::GetAdditionalSourceInf()
{
    std::vector<double> res(3);
    res[0] = Z1;
    res[1] = Z2;
    res[2] = ndZ;

    return res;
}

template <class PointType>
void EmitterDevice3d<PointType>::SetAdditionalSourceInf(std::vector<double> inf)
{
    Z1  = inf[0];
    Z2  = inf[1];
    ndZ = inf[2];
}

template <class PointType>
void EmitterDevice3d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData2d<PointType>>&                grid)
{
    this->boundaryList = in;
    this->particleSource->InitEmissionBoundary(boundaryIn, grid, parametersIn, error);
}

template <class PointType>
void EmitterDevice3d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData2dpolar<PointType>>&           grid)
{
    this->boundaryList = in;
    this->particleSource->InitEmissionBoundary(boundaryIn, grid, parametersIn, error);
}

template <class PointType>
void EmitterDevice3d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData3d<PointType>>&                grid)
{
    /*this->boundaryList = in;
    this->particleSources.clear();
    zArray.clear();
    this->particleSources.push_back( std::shared_ptr<ParticleSource2d<PointType>>(new
    ParticleSource2d<PointType>()));
    this->particleSources[0]->InitEmissionBoundary(boundaryIn, grid); double dz = (Z2 - Z1) / ndZ;
    this->particleSources.back()->SetZCoordinate(Z1 + dz / 2);

    double zCur = Z1;
    for (int i = 0; i < ndZ; i++)
    {
            zArray.push_back(zCur);

            if (i > 0)
            {
                    this->particleSources.push_back(this->particleSources[0]);
                    this->particleSources.back()->SetZCoordinate(zCur);
            }

            zCur = zCur + dz;

    };
    zArray.push_back(Z2);
    this->particleSources.push_back(this->particleSources[0]);
    this->particleSources.back()->SetZCoordinate(Z2);*/
}

template <class PointType>
void EmitterDevice3d<PointType>::GenerateSyncParticle(
    const std::shared_ptr<Particles3d<PointType>>& particlesData, PointType restMass){

}

template <class PointType>
void EmitterDevice3d<PointType>::SetFlowCurrent(double res)
{
    this->particleSources[0]->SetFlowCurrent(res);
}

template <class PointType>
void EmitterDevice3d<PointType>::SetBoundariesList(
    std::vector<int> in, std::vector<double> parametersIn, std::string& error,
    std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData3d<PointType>>&                grid){

}

template <class PointType>
template <class Archive>
void EmitterDevice3d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
    ar& nParticlesXY;
    ar& phiXY;
    ar& nParticlesZ;
    ar& Z1;
    ar& Z2;
    ar& ndZ;
    ar & this->particleSources;
    ar& zArray;
}
template <class PointType>
template <class Archive>
void EmitterDevice3d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
    ar& nParticlesXY;
    ar& phiXY;
    ar& nParticlesZ;
    ar& Z1;
    ar& Z2;
    ar& ndZ;
    ar & this->particleSources;
    ar& zArray;
}

template class EmitterDevice3d<float>;
template class EmitterDevice3d<double>;

template void
EmitterDevice3d<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);

template void EmitterDevice3d<double>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void
EmitterDevice3d<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);

template void EmitterDevice3d<float>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;