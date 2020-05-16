
#include "Particle.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "GridData.h"
#include "ParticleGridInterface.h"
#include "ParticleShape2dCIC.h"
#include "ParticleShape2dTSC.h"
#include "WcalculateVector.h"

template class ParticleGridInterface<float>;
template class ParticleGridInterface<double>;

template <class PointType>
void ParticleGridInterface<PointType>::Grid2Particles(
    const std::shared_ptr<Particles2dpolar<PointType>>& particles,
    particlesFields<PointType>&                         particlesTmp,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType (*W)[9], int i1, int i2,
    int step, int recalculate)
{
    int       indexLoc;
    PointType beta1, beta2, beta3;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Er(),
                                           particlesTmp.Get_Er()[i - i1]);
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ephi(),
                                           particlesTmp.Get_Ephi()[i - i1]);
        }
        break;
    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            if (indexLoc == -1)
                continue;
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Er(),
                                           particlesTmp.Get_Er()[i - i1]);
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ephi(),
                                           particlesTmp.Get_Ephi()[i - i1]);
        }
    }
};
template <class PointType>
void ParticleGridInterface<PointType>::Particles2Grid(
    const std::shared_ptr<Particles2dpolar<PointType>>& particles, PointType* rho,
    PointType (*W)[9], int i1, int i2){

};

template <class PointType>
void ParticleGridInterface<PointType>::Wcalculate(
    const std::shared_ptr<Particles2dpolar<PointType>>& particles, PointType (*W)[9], int i1,
    int i2, int thread, int flagStepEstimate)
{
    int indexLoc;
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = 0; i < particles->NParticles(); i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            for (int k  = 0; k < 9; k++)
                W[i][k] = 0;
            TSCcellArray->WcalculatePolar(indexLoc, particles->Get_r()[i], particles->Get_phi()[i],
                                          &W[i][0]);
        };
        break;
    case 0:

        searchIndexesLarge[thread].resize(i2 - i1);
        searchIndexesLarge1[thread].clear();
        for (int i = i1; i < i2; i++)
        {
            searchIndexesLarge[thread][i - i1] =
                templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            if (searchIndexesLarge[thread][i - i1] == -1)
            {
                searchIndexesLarge[thread][i - i1] = 1;
                searchIndexesLarge1[thread].push_back(i - i1);
            }
        }

        tmpData1[thread].resize(i2 - i1);
        tmpData2[thread].resize(i2 - i1);

        if (1 == flagStepEstimate)
        {
            particles->GammaCalc(gamma[thread]);
            CICcellArray->WcalculatePolar(
                &searchIndexesLarge[thread][0], &particles->Get_r()[0], &particles->Get_phi()[0], W,
                i1, i2, &particles->Get_pr()[0], &particles->Get_pphi()[0], &tmpData1[thread][0],
                &tmpData2[thread][0], &gamma[thread][0], particles->minStep);
        }
        else
            CICcellArray->WcalculatePolar(&searchIndexesLarge[thread][0], &particles->Get_r()[0],
                                          &particles->Get_phi()[0], W, i1, i2);

        for (int i = 0; i < searchIndexesLarge1[thread].size(); i++)
        {
            for (int k                               = 0; k < 9; k++)
                W[searchIndexesLarge1[thread][i]][k] = 0;
        };
        break;
    }
};
