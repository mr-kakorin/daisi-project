#include "ParticleGridInterface.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "FlagStringsSolver.h"
#include "GridData.h"
#include "Particle.h"
#include "ParticleShape2d.h"
#include "ParticleShape2dCIC.h"
#include "ParticleShape2dTSC.h"
#include "WcalculateVector.h"


template <class PointType>
long long ParticleGridInterface<PointType>::GetMemorySize()
{

    long long result = 0;
    result           = result + Dmath::vectorsize(CICArray);
    result           = result + Dmath::vectorsize(isBoundary);
    result           = result + Dmath::vectorsize(cellTypes);
    result           = result + Dmath::vectorsize(searchIndexes);
    result           = result + Dmath::vectorsize(searchIndexesLarge);
    result           = result + Dmath::vectorsize(searchIndexesLarge1);
    result           = result + Dmath::vectorsize(searchIndexesSpecial);
    result           = result + Dmath::vectorsize(tmpData1);
    result           = result + Dmath::vectorsize(tmpData2);
    result           = result + CICcellArray->GetMemorySize();
    return result;
}

template <class PointType>
void ParticleGridInterface<PointType>::init(
    int nflow, const std::shared_ptr<GridData3d<PointType>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer3d<PointType>> domainBoundary, int size, int numThreads,
    int blockSize){}

template <class PointType>
std::vector<unsigned int>& ParticleGridInterface<PointType>::CheckParticlesBoundaries(
    const std::shared_ptr<BoundaryConditions>&                          boundaryConditions,
    std::vector<unsigned int>&                                          result,
    const std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
    const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>&    conductorList,
    const std::shared_ptr<Particles3d<PointType>>&                      state1,
    const std::shared_ptr<Particles3d<PointType>>& state2, PointType dt, PointType charge,
    PointType mass, int i1, int i2, int thread)
{
    std::vector<unsigned int> t;
    return t;
}

template <class PointType>
void ParticleGridInterface<PointType>::InitEmCells(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    const std::shared_ptr<Particles3d<PointType>>& state1, double dH, int flag, int flagInit){}

template <class PointType>
void ParticleGridInterface<PointType>::Particles2GridPTI(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
    const std::shared_ptr<Particles3d<PointType>>& state1, PointType (*W1)[9],
    const std::shared_ptr<Particles3d<PointType>>& state2, PointType (*W2)[9], PointType* rho,
    PointType dt, int flagType, int i1, int i2, int thread){}

template <class PointType>
std::vector<double> ParticleGridInterface<PointType>::GetParameters()
{
    std::vector<double> result;
    result.push_back(double(epsilon));
    result.push_back(fabsorpFrac);
    result.push_back(priorityParticleShapeType);
    result.push_back(double(removeParticlesFlag));
    return result;
}

template <class PointType>
void ParticleGridInterface<PointType>::SetParameters(std::vector<double> in)
{
    priorityParticleShapeType = in[2];
    epsilon                   = in[0];
    fabsorpFrac                = in[1];
    removeParticlesFlag       = in[3];
}

template <class PointType>
ParticleGridInterface<PointType>::ParticleGridInterface()
{
    priorityParticleShapeType = 0;
    epsilon                   = 1e-5;
    fabsorpFrac                = 1e-3;
    removeParticlesFlag       = 1;
}

template <class PointType>
void ParticleGridInterface<PointType>::SearchStartCellsEmission(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
    std::vector<std::shared_ptr<Particles3d<PointType>>>&         particles){}

template <class PointType>
int ParticleGridInterface<PointType>::InCell(PointType x1, PointType x2)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray->InCell(x1, x2);
    case 0:
        return CICcellArray->InCell(x1, x2);
    }
    return -1;
}

template <class PointType>
int ParticleGridInterface<PointType>::InCellWithEps(PointType x1, PointType x2)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray->InCellWithEps(x1, x2);
    case 0:
        return CICcellArray->InCellWithEps(x1, x2);
    }
    return -1;
}

template <class PointType>
PointType ParticleGridInterface<PointType>::GetH1(int index)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray->GetH1(index);
    case 0:
        return CICcellArray->GetH1(index);
    }
}

template <class PointType>
PointType ParticleGridInterface<PointType>::GetH2(int index)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray->GetH2(index);
    case 0:
        return CICcellArray->GetH2(index);
    }
}

template <class PointType>
int ParticleGridInterface<PointType>::InCell(PointType x1, PointType x2,
                                             const std::vector<unsigned int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (TSCcellArray->InCell(index, x1, x2))
                    return searchIndexes[i];
            }
        }
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (CICcellArray->InCell(index, x1, x2))
                    return searchIndexes[i];
            }
        }
        break;
    }
    return -1;
}

template <class PointType>
int ParticleGridInterface<PointType>::InCellWithEps(PointType x1, PointType x2,
                                                    const std::vector<unsigned int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (TSCcellArray->InCellWithEps(index, x1, x2))
                    return searchIndexes[i];
            }
        }
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (CICcellArray->InCellWithEps(index, x1, x2))
                    return searchIndexes[i];
            }
        }
        break;
    }
    return -1;
}

template <class PointType>
int ParticleGridInterface<PointType>::InCellWithEps(PointType x1, PointType x2, PointType x3,
                                                    const std::vector<unsigned int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0, 0);
            if (index != -1)
            {
                if (CICcellArray->InCellWithEps(index, x1, x2, x3))
                    return searchIndexes[i];
            }
        }
        break;
    }
    return -1;
}

template <class PointType>
int ParticleGridInterface<PointType>::InCell(PointType x1, PointType x2, PointType x3,
                                             const std::vector<unsigned int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb->getElem(searchIndexes[i], 0, 0, 0);
            if (index != -1)
            {
                if (CICcellArray->InCellWithEps(index, x1, x2, x3))
                    return searchIndexes[i];
            }
        }
        break;
    }
    return -1;
}

/*int ParticleGridInterface <PointType>::InCellNew(PointType x1, PointType x2, const
std::vector<int>&
searchIndexes,const std::shared_ptr<GridData2daxs<DataContainer>* const gridData)
{
        switch (priorityParticleShapeType)
        {
        case 1:
                for (int i = 0; i < searchIndexes.size(); i++)
                {
                        int index = templNumb->getElem(searchIndexes[i], 0, 0);
                        if (index != -1)
                        {
                                if (cellArray[index]->InCell(x1, x2, gridData->Getr(),
gridData->Getz()))
                                        return searchIndexes[i];
                        }
                };
                break;
        case 0:
                for (int i = 0; i < searchIndexes.size(); i++)
                {
                        int index = templNumb->getElem(searchIndexes[i], 0, 0);
                        if (index != -1)
                        {
                                if (ParticleShapeCIC2dStatic <PointType>::InCell(x1, x2,
gridData->Getr(),
gridData->Getz(), index, CICArray[index])) return searchIndexes[i];
                        }
                };
                break;
        }
        return -1;
};*/

template <class PointType>
void ParticleGridInterface<PointType>::Grid2Particles(
    const std::shared_ptr<Particles2d<PointType>>&      particles,
    particlesFields<PointType>&                         particlesTmp,
    const std::shared_ptr<const GridData2d<PointType>>& gridData, PointType (*W)[9], int i1, int i2,
    int step, int recalculate)
{
    // gridData->Get_Ex()[0] = 1;
    if (i1 == i2)
        return;
    int       indexLoc;
    PointType beta1, beta2, beta3;

    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ex(),
                                           particlesTmp.Get_Ex()[i - i1]);
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ey(),
                                           particlesTmp.Get_Ey()[i - i1]);
        }
        break;
    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);

            if (indexLoc == -1)
            {
                particlesTmp.Get_Ex()[i - i1] = 0;
                particlesTmp.Get_Ey()[i - i1] = 0;
                continue;
            }

            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ex(),
                                           particlesTmp.Get_Ex()[i - i1]);
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ey(),
                                           particlesTmp.Get_Ey()[i - i1]);

            if (step % recalculate == 0)
            {
                CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_ExCol(),
                                               particles->Get_ExCol()[i]);
                CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_EyCol(),
                                               particles->Get_EyCol()[i]);
            }

            particlesTmp.Get_Ex()[i - i1] =
                particles->Get_ExCol()[i] + particlesTmp.Get_Ex()[i - i1];
            particlesTmp.Get_Ey()[i - i1] =
                particles->Get_EyCol()[i] + particlesTmp.Get_Ey()[i - i1];
        }
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Grid2Particles(
    const std::shared_ptr<Particles3d<PointType>>& particles,
    particlesFields<PointType>&                    particlesTmp,
    const std::shared_ptr<GridData3d<PointType>>& gridData, PointType (*W)[9], int i1, int i2,
    int step, int recalculate)
{
    if (i1 == i2)
        return;
    int       indexLoc;
    PointType beta1, beta2, beta3;

    switch (priorityParticleShapeType)
    {
    case 1:
        break;
    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);

            if (indexLoc == -1)
            {
                particlesTmp.Get_Ex()[i - i1] = 0;
                particlesTmp.Get_Ey()[i - i1] = 0;
                particlesTmp.Get_Ez()[i - i1] = 0;

                particlesTmp.Get_Bx()[i - i1] = 0;
                particlesTmp.Get_By()[i - i1] = 0;
                particlesTmp.Get_Bz()[i - i1] = 0;
                continue;
            }

            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->Get_Ex(),
                                             particlesTmp.Get_Ex()[i - i1]);
            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->Get_Ey(),
                                             particlesTmp.Get_Ey()[i - i1]);
            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->Get_Ez(),
                                             particlesTmp.Get_Ez()[i - i1]);

            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->GetBx(),
                                             particlesTmp.Get_Bx()[i - i1]);
            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->GetBy(),
                                             particlesTmp.Get_By()[i - i1]);
            CICcellArray->ValueInterpolate3d(indexLoc, &W[i - i1][0], gridData->GetBz(),
                                             particlesTmp.Get_Bz()[i - i1]);
        }
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Grid2Particles(
    const std::shared_ptr<Particles3dcil<PointType>>& particles,
    particlesFields<PointType>&                       particlesTmp,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType (*W)[9], int i1, int i2,
    int step, int recalculate)
{
    if (i1 == i2)
        return;
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
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ez(),
                                           particlesTmp.Get_Ez()[i - i1]);
            TSCcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->GetBphi(),
                                           particlesTmp.Get_Bphi()[i - i1]);
        }
        break;
    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            if (indexLoc == -1)
            {
                particlesTmp.Get_Er()[i - i1]   = 0;
                particlesTmp.Get_Ez()[i - i1]   = 0;
                particlesTmp.Get_Bphi()[i - i1] = 0;
                continue;
            }
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Er(),
                                           particlesTmp.Get_Er()[i - i1]);
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_Ez(),
                                           particlesTmp.Get_Ez()[i - i1]);
            CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->GetBphi(),
                                           particlesTmp.Get_Bphi()[i - i1]);

            if (step % recalculate == 0)
            {
                CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_ErCol(),
                                               particles->Get_ErCol()[i]);
                CICcellArray->ValueInterpolate(indexLoc, &W[i - i1][0], gridData->Get_EzCol(),
                                               particles->Get_EzCol()[i]);
            }

            particlesTmp.Get_Er()[i - i1] =
                particles->Get_ErCol()[i] + particlesTmp.Get_Er()[i - i1];
            particlesTmp.Get_Ez()[i - i1] =
                particles->Get_EzCol()[i] + particlesTmp.Get_Ez()[i - i1];

            //	particles->GetBetaComponents(beta1, beta2, beta3, i);
            //	particles->cellSize[i] = CICcellArray->GetSize(beta1, beta2, beta3, indexLoc);
        }
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Particles2Grid(
    const std::shared_ptr<Particles2d<PointType>>& particles, PointType* rho, PointType (*W)[9],
    int i1, int i2)
{
    if (i1 == i2)
        return;
    int indexLoc;

    for (int i = i1; i < i2; i++)
    {
        indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
        if (-1 == indexLoc)
            continue;
        CICcellArray->ChargeCalculate(indexLoc, &W[i - i1][0], particles->q[i], rho);
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Particles2Grid(
    const std::shared_ptr<Particles3dcil<PointType>>& particles, PointType* rho, PointType (*W)[9],
    int i1, int i2)
{

    if (i1 == i2)
        return;
    int indexLoc;

    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray->ChargeCalculate(indexLoc, &W[i - i1][0], particles->q[i], rho);
        }
        break;

    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            if (-1 == indexLoc)
                continue;
            CICcellArray->ChargeCalculate(indexLoc, &W[i - i1][0], particles->q[i], rho);
        }
        break;
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Particles2Grid(
    const std::shared_ptr<Particles3d<PointType>>& particles, PointType* rho, PointType (*W)[9],
    int i1, int i2)
{

    if (i1 == i2)
        return;
    int indexLoc;

    switch (priorityParticleShapeType)
    {
    case 1:
        break;

    case 0:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            if (-1 == indexLoc)
                continue;
            CICcellArray->ChargeCalculate3d(indexLoc, &W[i - i1][0], particles->q[i], rho);
        }
        break;
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::Particles2GridPTI(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
    const std::shared_ptr<particlesType>& state1, PointType (*W1)[9],
    const std::shared_ptr<particlesType>& state2, PointType (*W2)[9], PointType* rho, PointType dt,
    int flag, int i1, int i2, int thread)
{
    if (i1 == i2)
        return;
    int index1;
    int index2;

    PointType* Xar1 = state1->GetPointerToPosition1();
    PointType* Yar1 = state1->GetPointerToPosition2();

    PointType* Xar2 = state2->GetPointerToPosition1();
    PointType* Yar2 = state2->GetPointerToPosition2();

    PointType* CartArX1 = state1->GetPointerToCartesianX();
    PointType* CartArX2 = state1->GetPointerToCartesianY();

    PointType* CartArX12 = state2->GetPointerToCartesianX();
    PointType* CartArX22 = state2->GetPointerToCartesianY();
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = i1; i < i2; i++)
        {
            index1 = templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
            index2 = templNumb->getElem(state2->cellsNumbers[i], 0, 0);

            TSCcellArray->ChargeCalculate(Xar1[i], Yar1[i], index1, &W1[i - i1][0], Xar2[i],
                                          Yar2[i], index2, &W2[i - i1][0], 0, 0, state2->q[i], dt,
                                          rho);
        }
        break;

    case 0:
        /*for (int i = i1; i < i2; i++)
        {
                index1 = templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
                index2 = templNumb->getElem(state2->cellsNumbers[i], 0, 0);

                if (-1 == index1)
                        continue;

                if (-1 == index2)
                        index2 = index1;
        }*/

        searchIndexesLarge[thread].resize(i2 - i1);
        searchIndexesLarge1[thread].resize(i2 - i1);

        for (int i = i1; i < i2; i++)
        {
            searchIndexesLarge[thread][i - i1] =
                templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
            //	if (searchIndexesLarge[thread][i - i1] == -1)
            //		searchIndexesLarge[thread][i - i1] = 1;

            searchIndexesLarge1[thread][i - i1] = templNumb->getElem(state2->cellsNumbers[i], 0, 0);
            if (searchIndexesLarge1[thread][i - i1] == -1)
                searchIndexesLarge1[thread][i - i1] = searchIndexesLarge[thread][i - i1];
        }

        unsigned int* p;
        CICcellArray->ChargeCalculate(
            p, i1, i2, nearCathodeVolumes, emType, Icoef, &state2->Get_currentCoef()[0],
            &state2->flagEmitted[0], emissionCells, &state2->Get_startCellNumbers()[0], Xar1, Yar1,
            &searchIndexesLarge[thread][0], W1, Xar2, Yar2, &searchIndexesLarge1[thread][0], W2,
            &state2->q[0], dt, rho, CartArX1, CartArX2, flag, CartArX12, CartArX22);

        break;
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Particles2GridPTI(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
    const std::shared_ptr<Particles2dpolar<PointType>>& state1, PointType (*W1)[9],
    const std::shared_ptr<Particles2dpolar<PointType>>& state2, PointType (*W2)[9], PointType* rho,
    PointType dt, int flag, int i1, int i2, int thread)
{
    if (i1 == i2)
        return;
    int index1;
    int index2;

    PointType* Xar1 = state1->GetPointerToPosition1();
    PointType* Yar1 = state1->GetPointerToPosition2();

    PointType* Xar2 = state2->GetPointerToPosition1();
    PointType* Yar2 = state2->GetPointerToPosition2();

    PointType* CartArX1 = state1->GetPointerToCartesianX();
    PointType* CartArX2 = state1->GetPointerToCartesianY();

    PointType* CartArX12 = state2->GetPointerToCartesianX();
    PointType* CartArX22 = state2->GetPointerToCartesianY();
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = i1; i < i2; i++)
        {
            index1 = templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
            index2 = templNumb->getElem(state2->cellsNumbers[i], 0, 0);

            TSCcellArray->ChargeCalculate(Xar1[i], Yar1[i], index1, &W1[i - i1][0], Xar2[i],
                                          Yar2[i], index2, &W2[i - i1][0], 0, 0, state2->q[i], dt,
                                          rho);
        }
        break;

    case 0:
        /*for (int i = i1; i < i2; i++)
        {
        index1 = templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
        index2 = templNumb->getElem(state2->cellsNumbers[i], 0, 0);

        if (-1 == index1)
        continue;

        if (-1 == index2)
        index2 = index1;
        }*/

        searchIndexesLarge[thread].resize(i2 - i1);
        searchIndexesLarge1[thread].resize(i2 - i1);

        for (int i = i1; i < i2; i++)
        {
            searchIndexesLarge[thread][i - i1] =
                templNumb->getElem(state1->cellsNumbers[i - i1], 0, 0);
            if (searchIndexesLarge[thread][i - i1] == -1)
                searchIndexesLarge[thread][i - i1] = 1;

            searchIndexesLarge1[thread][i - i1] = templNumb->getElem(state2->cellsNumbers[i], 0, 0);
            if (searchIndexesLarge1[thread][i - i1] == -1)
                searchIndexesLarge1[thread][i - i1] = searchIndexesLarge[thread][i - i1];
        }

        CICcellArray->ChargeCalculate(
            &state2->Get_isPeriodical()[0], i1, i2, nearCathodeVolumes, emType, Icoef,
            &state2->Get_currentCoef()[0], &state2->flagEmitted[0], emissionCells,
            &state2->Get_startCellNumbers()[0], Xar1, Yar1, &searchIndexesLarge[thread][0], W1,
            Xar2, Yar2, &searchIndexesLarge1[thread][0], W2, &state2->q[0], dt, rho, CartArX1,
            CartArX2, flag, CartArX12, CartArX22);

        break;
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::InitEmCells(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    const std::shared_ptr<particlesType>& state1, double dH, int flag, int flagInit)
{

    PointType* Xar = state1->GetPointerToCartesianX();
    PointType* Yar = state1->GetPointerToCartesianY();

    int nrow = templNumb->GetNrow();

    std::vector<unsigned int> searchIndexesL(9);

    for (int i = 0; i < state1->NParticles(); i++)
    {
        state1->Get_startCellNumbers()[i] = -1;
        PointType px                      = state1->GetCartesianPX(i);
        PointType py                      = state1->GetCartesianPY(i);
        PointType p                       = sqrt(px * px + py * py);
        PointType x1                      = Xar[i];
        PointType y1                      = Yar[i];
        PointType xtmp;
        PointType ytmp;

        x1 = x1 + dH * 0.0002 * px / p;
        y1 = y1 + dH * 0.0002 * py / p;

        for (int j = 0; j < 11; j++)
        {
            int base = state1->cellsNumbers[i];

            searchIndexesL[0] = base;
            searchIndexesL[1] = base - 1;
            searchIndexesL[2] = base + 1;

            searchIndexesL[3] = base + nrow;
            searchIndexesL[4] = base - 1 + nrow;
            searchIndexesL[5] = base + 1 + nrow;

            searchIndexesL[6] = base - nrow;
            searchIndexesL[7] = base - 1 - nrow;
            searchIndexesL[8] = base + 1 - nrow;

            if (flag == 3)
            {
                Dmath::Cartesian2Polar(x1, y1, xtmp, ytmp);
            }
            else
            {
                xtmp = x1;
                ytmp = y1;
            }

            int index = InCellWithEps(xtmp, ytmp, searchIndexesL);

            if (index != -1)
            {
                int index1 = templNumb->getElem(index, 0, 0);

                int n = CICcellArray->InitEmCells(nearCathodeVolumes, x1, y1, index1, flagInit);
                if (state1->Get_startCellNumbers()[i] == -1)
                    state1->Get_startCellNumbers()[i] = n;
            }
            x1 = x1 + dH * 0.1 * px / p;
            y1 = y1 + dH * 0.1 * py / p;
        }
        if (state1->Get_startCellNumbers()[i] == -1 && i > 0)
            state1->Get_startCellNumbers()[i] = state1->Get_startCellNumbers()[i - 1];
    }

    for (int i = 0; i < state1->NParticles(); i++)
    {
        if (state1->Get_startCellNumbers()[i] == -1)
        {
            int j = 1;
            while (1)
            {
                if (j + i >= state1->Get_startCellNumbers().size())
                    break;
                if (state1->Get_startCellNumbers()[i + j] != -1)
                {
                    state1->Get_startCellNumbers()[i] = state1->Get_startCellNumbers()[i + j];
                    break;
                }
                j++;
            }
        }
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Charge2Density(PointType* rho, std::vector<int>& nonZeros)
{
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int j = 0; j < TSCcellArray->cellVolume.size(); j++)
        {
            if (rho[j])
                nonZeros.push_back(j);

            rho[j] = rho[j] / TSCcellArray->CellVolume(j);
        }
        break;
    case 0:
        for (int j = 0; j < CICcellArray->cellVolume.size(); j++)
        {
            // rho[j] = rho[j] / CICcellArray->CellVolume(j);
            if (rho[j])
            {
                nonZeros.push_back(j);
                rho[j] = rho[j] / CICcellArray->CellVolume(j);
                // rho[j] = -0.05;
            }
            //	if (std::isnan(rho[j]) || std::isinf(rho[j]))
            //		rho[j] = 0;
        }
        break;
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Wcalculate(
    const std::shared_ptr<Particles3dcil<PointType>>& particles, PointType (*W)[9], int i1, int i2,
    int thread, int flagStepEstimate)
{
    if (i2 == i1)
        return;

    int indexLoc;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            for (int k       = 0; k < 9; k++)
                W[i - i1][k] = 0;
            TSCcellArray->Wcalculate(indexLoc, particles->Get_r()[i], particles->Get_z()[i],
                                     &W[i - i1][0]);
        }
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
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_r()[0],
                                     &particles->Get_z()[0], W, i1, i2, &particles->Get_pr()[0],
                                     &particles->Get_pz()[0], &tmpData1[thread][0],
                                     &tmpData2[thread][0], &gamma[thread][0], particles->minStep);
        }
        else
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_r()[0],
                                     &particles->Get_z()[0], W, i1, i2);

        for (int i = 0; i < searchIndexesLarge1[thread].size(); i++)
        {
            for (int k                               = 0; k < 9; k++)
                W[searchIndexesLarge1[thread][i]][k] = 0;
        }

        break;
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Wcalculate(
    const std::shared_ptr<Particles2d<PointType>>& particles, PointType (*W)[9], int i1, int i2,
    int thread, int flagStepEstimate)
{
    if (i1 == i2)
        return;
    int indexLoc;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = i1; i < i2; i++)
        {
            indexLoc = templNumb->getElem(particles->cellsNumbers[i], 0, 0);
            for (int k       = 0; k < 9; k++)
                W[i - i1][k] = 0;
            TSCcellArray->Wcalculate(indexLoc, particles->Get_x()[i], particles->Get_y()[i],
                                     &W[i][0]);
        }
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
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_x()[0],
                                     &particles->Get_y()[0], W, i1, i2, &particles->Get_px()[0],
                                     &particles->Get_py()[0], &tmpData1[thread][0],
                                     &tmpData2[thread][0], &gamma[thread][0], particles->minStep);
        }
        else
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_x()[0],
                                     &particles->Get_y()[0], W, i1, i2);

        for (int i = 0; i < searchIndexesLarge1[thread].size(); i++)
        {
            for (int k                               = 0; k < 9; k++)
                W[searchIndexesLarge1[thread][i]][k] = 0;
        }

        break;
    }
}

template <class PointType>
void ParticleGridInterface<PointType>::Wcalculate(
    const std::shared_ptr<Particles3d<PointType>>& particles, PointType (*W)[9], int i1, int i2,
    int thread, int flagStepEstimate)
{
    if (i1 == i2)
        return;
    switch (priorityParticleShapeType)
    {
    case 1:
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
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_x()[0],
                                     &particles->Get_y()[0], &particles->Get_z()[0], W, i1, i2,
                                     &particles->Get_px()[0], &particles->Get_py()[0],
                                     &particles->Get_pz()[0], &tmpData1[thread][0],
                                     &tmpData2[thread][0], &gamma[thread][0], particles->minStep);
        }
        else
            CICcellArray->Wcalculate(&searchIndexesLarge[thread][0], &particles->Get_x()[0],
                                     &particles->Get_y()[0], &particles->Get_z()[0], W, i1, i2);

        for (int i = 0; i < searchIndexesLarge1[thread].size(); i++)
        {
            for (int k                               = 0; k < 9; k++)
                W[searchIndexesLarge1[thread][i]][k] = 0;
        }
        break;
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::SearchStartCells(
    int flow, const std::shared_ptr<particlesType>& particles)
{
    startCellNumbers[flow].clear();
    std::vector<unsigned int> searchIndexesLoc;

    searchIndexesLoc.resize(templNumb->ArSize());
    for (int k              = 0; k < templNumb->ArSize(); k++)
        searchIndexesLoc[k] = k;

    particles->cellsNumbers.clear();
    particles->cellsNumbers.resize(particles->NParticles());

    startCellNumbers[flow].clear();
    startCellNumbers[flow].resize(particles->NParticles());

    for (int i = 0; i < particles->NParticles(); i++)
    {
        particles->cellsNumbers[i] =
            InCellWithEps(particles->GetPointerToPosition1()[i],
                          particles->GetPointerToPosition2()[i], searchIndexesLoc);

        if (particles->flagEmitted[i] == 0)
        {
            startCellNumbers[flow][i] =
                InCellWithEps(particles->GetPointerToPosition1()[i],
                              particles->GetPointerToPosition2()[i], searchIndexesLoc);
            particles->cellsNumbers[i] = startCellNumbers[flow][i];
        }
    }

    /*std::vector<int>  searchIndexesParticles (startCellNumbers[flow].size()*9);

    int nrow = templNumb->GetNrow();
    int base;

    for (int j = 0; j < startCellNumbers[flow].size(); j++)
    {
            base = startCellNumbers[flow][j];
            searchIndexesParticles[9 * j + 0] = base;
            searchIndexesParticles[9 * j + 1] = base - 1;
            searchIndexesParticles[9 * j + 2] = base + 1;

            searchIndexesParticles[9 * j + 3] = base + nrow;
            searchIndexesParticles[9 * j + 4] = base - 1 + nrow;
            searchIndexesParticles[9 * j + 5] = base + 1 + nrow;

            searchIndexesParticles[9 * j + 6] = base - nrow;
            searchIndexesParticles[9 * j + 7] = base - 1 - nrow;
            searchIndexesParticles[9 * j + 8] = base + 1 - nrow;
    }
    searchIndexes.resize(9);
    */
}

template <class PointType>
void ParticleGridInterface<PointType>::SearchStartCells(
    int flow, const std::shared_ptr<Particles3d<PointType>>& particles)
{
    startCellNumbers[flow].clear();
    std::vector<unsigned int> searchIndexesLoc;

    searchIndexesLoc.resize(templNumb->ArSize());
    for (int k              = 0; k < templNumb->ArSize(); k++)
        searchIndexesLoc[k] = k;

    particles->cellsNumbers.clear();
    particles->cellsNumbers.resize(particles->NParticles());

    startCellNumbers[flow].clear();
    startCellNumbers[flow].resize(particles->NParticles());

    for (int i = 0; i < particles->NParticles(); i++)
    {
        particles->cellsNumbers[i] = InCellWithEps(
            particles->GetPointerToPosition1()[i], particles->GetPointerToPosition2()[i],
            particles->GetPointerToPosition3()[i], searchIndexesLoc);

        if (particles->flagEmitted[i] == 0)
        {
            startCellNumbers[flow][i] = InCellWithEps(
                particles->GetPointerToPosition1()[i], particles->GetPointerToPosition2()[i],
                particles->GetPointerToPosition3()[i], searchIndexesLoc);
            particles->cellsNumbers[i] = startCellNumbers[flow][i];
        }
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::SearchStartCellsEmission(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
    std::vector<std::shared_ptr<particlesType>>&                  particles)
{
    for (int i = 0; i < particles.size(); i++)
        for (int j = 0; j < particles[i]->Get_startCellNumbers().size(); j++)
            particles[i]->Get_currentCoef()[j] = -1;

    for (int i = 0; i < particles.size(); i++)
    {
        for (int j = 0; j < particles[i]->Get_startCellNumbers().size(); j++)
        {
            PointType val = particles[i]->Get_startCellNumbers()[j];
            PointType        currentFrom_dl = 0;
            std::vector<int> indexes1;
            std::vector<int> indexes2;
            for (int i1 = 0; i1 < particles.size(); i1++)
            {
                for (int j1 = 0; j1 < particles[i1]->Get_startCellNumbers().size(); j1++)
                {
                    if (val == particles[i1]->Get_startCellNumbers()[j1] &&
                        particles[i1]->Get_currentCoef()[j1] == -1)
                    {
                        currentFrom_dl = currentFrom_dl + particles[i1]->q[j1];
                        indexes1.push_back(i1);
                        indexes2.push_back(j1);
                    }
                }
            }
            for (int i1 = 0; i1 < indexes1.size(); i1++)
            {
                particles[indexes1[i1]]->Get_currentCoef()[indexes2[i1]] =
                    particles[indexes1[i1]]->q[indexes2[i1]] / std::abs(currentFrom_dl);
                if (std::abs(currentFrom_dl) < 1e-17)
                    particles[indexes1[i1]]->Get_currentCoef()[indexes2[i1]] = 0;
            }
        }
    }
}

template <class PointType>
template <class particlesType>
std::vector<unsigned int>
ParticleGridInterface<PointType>::InCell(const std::shared_ptr<particlesType>& particles,
                                         std::vector<unsigned int> EmptyPlaces, int i1, int i2,
                                         int thread, int flow)
{

    std::vector<unsigned int> remove;

    if (i1 == i2)
        return remove;

    int nrow = templNumb->GetNrow();
    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();
    int base;

    for (int i = i1; i < i2; i++)
    {
        if (particles->cellsNumbers[i] == -1)
            continue;

        int flagContinue = 0;
        /*	for (int j = 0; j < EmptyPlaces.size(); j++)
                {
                        if (EmptyPlaces[j] == i)
                        {
                                particles->cellsNumbers[i] = -1;
                                flagContinue = 1;
                                break;
                        }
                }*/

        if (flagContinue)
            continue;

        if (particles->flagEmitted[i] < 2 || !particles->cellsNumbers[i])
        {
            particles->cellsNumbers[i] = InCellWithEps(Xar[i], Yar[i], startCellNumbers[flow]);

            if (-1 == particles->cellsNumbers[i])
            {
                searchIndexesLarge[thread].resize(templNumb->ArSize());

                for (int k                        = 0; k < templNumb->ArSize(); k++)
                    searchIndexesLarge[thread][k] = k;

                particles->cellsNumbers[i] =
                    InCellWithEps(Xar[i], Yar[i], searchIndexesLarge[thread]);
            }
        }
        else
        {
            base                     = particles->cellsNumbers[i];
            searchIndexes[thread][0] = base;
            searchIndexes[thread][1] = base - 1;
            searchIndexes[thread][2] = base + 1;

            searchIndexes[thread][3] = base + nrow;
            searchIndexes[thread][4] = base - 1 + nrow;
            searchIndexes[thread][5] = base + 1 + nrow;

            searchIndexes[thread][6] = base - nrow;
            searchIndexes[thread][7] = base - 1 - nrow;
            searchIndexes[thread][8] = base + 1 - nrow;

            int res = InCell(Xar[i], Yar[i], searchIndexes[thread]);

            if (removeParticlesFlag == 1)
            {
                particles->cellsNumbers[i] = res;
                continue;
            }

            if (-1 == res)
            {

                volatile double tt     = Xar[i];
                volatile int    k      = 2;
                int             indexS = 0;
                while (-1 == res)
                {
                    indexS = 0;
                    searchIndexesLarge[thread].resize(pow(double(2 * k + 1), 2));

                    for (int s0 = -k; s0 <= k; s0++)
                    {
                        for (int s = -k; s <= k; s++)
                        {
                            int ind = base - s + s0 * nrow;
                            if (ind < 0)
                                ind = templNumb->ArSize() + ind;
                            searchIndexesLarge[thread].push_back(ind);
                            indexS++;
                        }
                    }
                    res = InCell(Xar[i], Yar[i], searchIndexesLarge[thread]);
                    k++;
                    if (k > 20 && -1 == res)
                        break;
                }

                if (-1 == res)
                {
                    searchIndexesLarge[thread].resize(templNumb->ArSize());

                    for (int k                        = 0; k < templNumb->ArSize(); k++)
                        searchIndexesLarge[thread][k] = k;

                    res = InCell(Xar[i], Yar[i], searchIndexesLarge[thread]);
                }
            }
            particles->cellsNumbers[i] = res;
        }
    }
    //	particles->removeParticle(remove);
    return remove;
}

template <class PointType>
template <class particlesType>
std::vector<unsigned int>
ParticleGridInterface<PointType>::InCellWithEps(const std::shared_ptr<particlesType>& particles,
                                                int i1, int i2, int thread, int flow)
{
    std::vector<unsigned int> remove;

    if (i1 == i2)
        return remove;

    int nrow = templNumb->GetNrow();

    int base;
    int k = 0;

    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();

    for (int i = i1; i < i2; i++)
    {

        if (particles->cellsNumbers[i] == -1)
            continue;

        if (particles->cellsNumbers[i] == 0)
        {
            k++;
            particles->cellsNumbers[i] = InCellWithEps(Xar[i], Yar[i], startCellNumbers[flow]);
        }
        else
        {

            base                     = particles->cellsNumbers[i];
            searchIndexes[thread][0] = base;
            searchIndexes[thread][1] = base - 1;
            searchIndexes[thread][2] = base + 1;

            searchIndexes[thread][3] = base + nrow;
            searchIndexes[thread][4] = base - 1 + nrow;
            searchIndexes[thread][5] = base + 1 + nrow;

            searchIndexes[thread][6] = base - nrow;
            searchIndexes[thread][7] = base - 1 - nrow;
            searchIndexes[thread][8] = base + 1 - nrow;

            particles->cellsNumbers[i] = InCellWithEps(Xar[i], Yar[i], searchIndexes[thread]);
        }

        if (particles->cellsNumbers[i] == -1)
        {
            double r = sqrt(Xar[i] * Xar[i] + Yar[i] * Yar[i]);
            remove.push_back(i);
            //	particles->cellsNumbers[i] = base;
        }
    }
    return remove;
}

template <class PointType>
std::vector<unsigned int> ParticleGridInterface<PointType>::InCellWithEps(
    const std::shared_ptr<Particles3d<PointType>>& particles, int i1, int i2, int thread, int flow)
{
    std::vector<unsigned int> remove;

    if (i1 == i2)
        return remove;

    int nrow = templNumb->GetNrow();
    int ncol = templNumb->GetNcol();

    int base;
    int k = 0;

    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();
    PointType* Zar = particles->GetPointerToPosition3();

    for (int i = i1; i < i2; i++)
    {

        if (particles->cellsNumbers[i] == -1)
            continue;

        if (particles->cellsNumbers[i] == 0)
        {
            k++;
            particles->cellsNumbers[i] =
                InCellWithEps(Xar[i], Yar[i], Zar[i], startCellNumbers[flow]);
        }
        else
        {
            base = particles->cellsNumbers[i];
            for (int iz = -1; iz <= 1; iz++)
            {
                searchIndexes[thread][0] = base + iz * nrow * ncol;
                searchIndexes[thread][1] = base - 1 + iz * nrow * ncol;
                searchIndexes[thread][2] = base + 1 + iz * nrow * ncol;

                searchIndexes[thread][3] = base + nrow + iz * nrow * ncol;
                searchIndexes[thread][4] = base - 1 + nrow + iz * nrow * ncol;
                searchIndexes[thread][5] = base + 1 + nrow + iz * nrow * ncol;

                searchIndexes[thread][6] = base - nrow + iz * nrow * ncol;
                searchIndexes[thread][7] = base - 1 - nrow + iz * nrow * ncol;
                searchIndexes[thread][8] = base + 1 - nrow + iz * nrow * ncol;

                particles->cellsNumbers[i] =
                    InCellWithEps(Xar[i], Yar[i], Zar[i], searchIndexes[thread]);

                if (particles->cellsNumbers[i] != -1)
                    break;
            }
        }

        if (particles->cellsNumbers[i] == -1)
        {
            remove.push_back(i);
            //	particles->cellsNumbers[i] = base;
        }
    }
    return remove;
}

template <class PointType>
std::vector<unsigned int>
ParticleGridInterface<PointType>::InCell(const std::shared_ptr<Particles3d<PointType>>& particles,
                                         std::vector<unsigned int> EmptyPlaces, int i1, int i2,
                                         int thread, int flow)
{
    std::vector<unsigned int> remove;

    if (i1 == i2)
        return remove;

    int nrow = templNumb->GetNrow();
    int ncol = templNumb->GetNcol();

    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();
    PointType* Zar = particles->GetPointerToPosition3();

    int base;

    for (int i = i1; i < i2; i++)
    {

        if (particles->cellsNumbers[i] == -1)
            continue;

        if (particles->flagEmitted[i] < 2)
        {
            particles->cellsNumbers[i] =
                InCellWithEps(Xar[i], Yar[i], Zar[i], startCellNumbers[flow]);
        }
        else
        {

            int res;
            base = particles->cellsNumbers[i];

            for (int iz = -1; iz <= 1; iz++)
            {
                searchIndexes[thread][0] = base + iz * nrow * ncol;
                searchIndexes[thread][1] = base - 1 + iz * nrow * ncol;
                searchIndexes[thread][2] = base + 1 + iz * nrow * ncol;

                searchIndexes[thread][3] = base + nrow + iz * nrow * ncol;
                searchIndexes[thread][4] = base - 1 + nrow + iz * nrow * ncol;
                searchIndexes[thread][5] = base + 1 + nrow + iz * nrow * ncol;

                searchIndexes[thread][6] = base - nrow + iz * nrow * ncol;
                searchIndexes[thread][7] = base - 1 - nrow + iz * nrow * ncol;
                searchIndexes[thread][8] = base + 1 - nrow + iz * nrow * ncol;

                res = InCellWithEps(Xar[i], Yar[i], Zar[i], searchIndexes[thread]);
                if (res != -1)
                    break;
            }

            if (removeParticlesFlag == 1)
            {
                particles->cellsNumbers[i] = res;
                continue;
            }

            if (-1 == res)
            {
                int k      = 2;
                int indexS = 0;
                for (int iz = -1; iz <= 1; iz++)
                {
                    while (-1 == res)
                    {

                        indexS = 0;
                        searchIndexesLarge[thread].resize(pow(double(2 * k + 1), 2));
                        searchIndexesLarge[thread].clear();

                        for (int s0 = -k; s0 <= k; s0++)
                        {
                            for (int s = -k; s <= k; s++)
                            {
                                int ind = base - s + s0 * nrow + iz * nrow * ncol;
                                if (ind < 0)
                                    ind                            = templNumb->ArSize() + ind;
                                searchIndexesLarge[thread][indexS] = ind;
                                indexS++;
                            }
                        }
                        res = InCell(Xar[i], Yar[i], Zar[i], searchIndexesLarge[thread]);
                        k++;
                        if (k > 20 && -1 == res)
                            break;
                    }
                    if (res != -1)
                        break;
                }

                if (-1 == res)
                {
                    searchIndexesLarge[thread].resize(templNumb->ArSize());

                    for (int k                        = 0; k < templNumb->ArSize(); k++)
                        searchIndexesLarge[thread][k] = k;

                    res = InCell(Xar[i], Yar[i], Zar[i], searchIndexesLarge[thread]);
                }
            }
            particles->cellsNumbers[i] = res;
        }

        if (particles->cellsNumbers[i] == 0)
        {
            int tt = 0;
        }

        if (particles->cellsNumbers[i] == -1)
        {
            remove.push_back(i);
            //	particles->cellsNumbers[i] = base;
        }
    }
    //	particles->removeParticle(remove);
    return remove;
}


template <class PointType>
template <class PointType1>
void ParticleGridInterface<PointType>::axsPolar(
    Particles3dcil<PointType1>* particles, const std::shared_ptr<GridData2dpolar<double>>& gridData,
    std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread){
    /*if (i1 == i2)
            return;

    std::vector<unsigned int> remove;
    int nrow = templNumb->GetNrow();

    int k = 0;

    int base;
    PointType W[4];

    PointType phi;

    double dphi = 0.058177641733144;
    double border1 = 1.5 * PI() - dphi / 2;
    double border2 = 1.5 * PI() + dphi / 2;

    int res;

    for (int i = i1; i < i2; i++)
    {
            if (particles->cellsNumbers[i] == -1)
                    continue;


            if (particles->Get_r()[i]>0.11 || particles->Get_r()[i]<0.09)
            {
    //		particles->Ephi[i] = 0;
                    continue;
            };

            phi = particles->Get_phi()[i] + dphi / 2;

            if (phi < border1)
            {
                    int c = (border1 - phi) / dphi + 1;
                    phi = phi + dphi * c;
            }


            if (phi > border2)
            {
                    int c = (phi - border2) / dphi + 1;
                    phi = phi - dphi * c;
            }


            if (phi < 1.5 * PI())
                    phi = 1.5 * PI() + (1.5 * PI() - phi);


            if (particles->cellsNumbersAdd[i] == 0)
                    particles->cellsNumbersAdd[i] = InCell(particles->Get_r()[i], phi,
    searchIndexesSpecial[thread]);

            else
            {
                    int base = particles->cellsNumbersAdd[i];

                    searchIndexes[thread][0] = base;
                    searchIndexes[thread][1] = base - 1;
                    searchIndexes[thread][2] = base + 1;

                    searchIndexes[thread][3] = base + nrow;
                    searchIndexes[thread][4] = base - 1 + nrow;
                    searchIndexes[thread][5] = base + 1 + nrow;

                    searchIndexes[thread][6] = base - nrow;
                    searchIndexes[thread][7] = base - 1 - nrow;
                    searchIndexes[thread][8] = base + 1 - nrow;

                    res = InCell(particles->Get_r()[i], phi, searchIndexes[thread]);


                    int k = 5;
                    res = -1;
                    while (-1 == res)
                    {
                            int indexS = 0;
                            searchIndexesLarge[thread].resize(pow(double(2 * k + 1), 2));
                            searchIndexesLarge[thread].clear();

                            for (int s0 = -1; s0 <= 1; s0++)
                            {
                                    for (int s = -k; s <= k; s++)
                                    {
                                            int ind = base - s + s0 * nrow;
                                            if (ind < 0)
                                                    ind = templNumb->ArSize() + ind;
                                            searchIndexesLarge[thread].push_back(ind);
                                            indexS++;
                                    }
                            }
                            res = InCell(particles->Get_r()[i], phi, searchIndexesLarge[thread]);
                            k++;
                            if (k > 6 && -1 == res)
                                    break;
                    };
                    particles->cellsNumbersAdd[i] = res;
            }



            if (particles->cellsNumbersAdd[i] == -1)
            {
            //	particles->Ephi[i] = 0;
                    continue;

            };


            int basePoint = templNumb->getElem(particles->cellsNumbersAdd[i], 0, 0);


            CICcellArray->WcalculatePolar(basePoint, particles->Get_r()[i], phi, W);

            PointType E;

            CICcellArray->ValueInterpolate(basePoint, W, gridData->Get_Ephi(), E);

    //	particles->Ephi[i] = E;
    };*/
}

/*template <class PointType>
std::vector<unsigned int> ParticleGridInterface <PointType>::InCell(const
std::shared_ptr<Particles3d<PointType>>& particles,
std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread)
{
        std::vector<unsigned int> remove;
        return remove;

};


template <class PointType>
std::vector<unsigned int> ParticleGridInterface <PointType>::InCellWithEps(const
std::shared_ptr<Particles3d<PointType>>&
particles, std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread)
{
        std::vector<unsigned int> remove;
        return remove;

};*/


template class ParticleGridInterface<float>;
template class ParticleGridInterface<double>;

template void ParticleGridInterface<float>::Particles2GridPTI<Particles3dcil<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
        arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles3dcil<float>>& state1, float (*W1)[9],
        const std::shared_ptr<Particles3dcil<float>>& state2, float (*W2)[9], float* rho, float dt,
        int flag, int i1, int i2, int thread);

template void ParticleGridInterface<float>::Particles2GridPTI<Particles2d<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
        arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles2d<float>>& state1, float (*W1)[9],
        const std::shared_ptr<Particles2d<float>>& state2, float (*W2)[9], float* rho, float dt,
        int flag, int i1, int i2, int thread);

template void ParticleGridInterface<double>::Particles2GridPTI<Particles3dcil<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
        arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles3dcil<double>>& state1, double (*W1)[9],
        const std::shared_ptr<Particles3dcil<double>>& state2, double (*W2)[9], double* rho, double dt,
        int flag, int i1, int i2, int thread);

template void ParticleGridInterface<double>::Particles2GridPTI<Particles2d<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
        arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles2d<double>>& state1, double (*W1)[9],
        const std::shared_ptr<Particles2d<double>>& state2, double (*W2)[9], double* rho, double dt,
        int flag, int i1, int i2, int thread);

template void ParticleGridInterface<float>::InitEmCells<Particles3dcil<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles3dcil<float>>& state1, double dH, int flag, int flagInit);
template void ParticleGridInterface<float>::InitEmCells<Particles2d<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles2d<float>>& state1, double dH, int flag, int flagInit);
template void ParticleGridInterface<float>::InitEmCells<Particles2dpolar<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles2dpolar<float>>& state1, double dH, int flag, int flagInit);

template void ParticleGridInterface<double>::InitEmCells<Particles3dcil<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles3dcil<double>>& state1, double dH, int flag, int flagInit);
template void ParticleGridInterface<double>::InitEmCells<Particles2d<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles2d<double>>& state1, double dH, int flag, int flagInit);
template void ParticleGridInterface<double>::InitEmCells<Particles2dpolar<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
        const std::shared_ptr<Particles2dpolar<double>>& state1, double dH, int flag, int flagInit);

template void ParticleGridInterface<float>::SearchStartCells<Particles3dcil<float>>(
        int flow, const std::shared_ptr<Particles3dcil<float>>& particles);

template void ParticleGridInterface<double>::SearchStartCells<Particles3dcil<double>>(
        int flow, const std::shared_ptr<Particles3dcil<double>>& particles);

template void ParticleGridInterface<float>::SearchStartCells<Particles2d<float>>(
        int flow, const std::shared_ptr<Particles2d<float>>& particles);

template void ParticleGridInterface<double>::SearchStartCells<Particles2d<double>>(
        int flow, const std::shared_ptr<Particles2d<double>>& particles);

template void ParticleGridInterface<float>::SearchStartCells<Particles2dpolar<float>>(
        int flow, const std::shared_ptr<Particles2dpolar<float>>& particles);

template void ParticleGridInterface<double>::SearchStartCells<Particles2dpolar<double>>(
        int flow, const std::shared_ptr<Particles2dpolar<double>>& particles);


template void ParticleGridInterface<float>::SearchStartCellsEmission<Particles3dcil<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles3dcil<float>>>&      particles);

template void ParticleGridInterface<double>::SearchStartCellsEmission<Particles3dcil<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles3dcil<double>>>&      particles);

template void ParticleGridInterface<float>::SearchStartCellsEmission<Particles2d<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles2d<float>>>&         particles);

template void ParticleGridInterface<double>::SearchStartCellsEmission<Particles2d<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles2d<double>>>&         particles);

template void ParticleGridInterface<float>::SearchStartCellsEmission<Particles2dpolar<float>>(
        const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles2dpolar<float>>>&    particles);

template void ParticleGridInterface<double>::SearchStartCellsEmission<Particles2dpolar<double>>(
        const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles2dpolar<double>>>&    particles);

template std::vector<unsigned int> ParticleGridInterface<float>::InCell<Particles3dcil<float>>(
        const std::shared_ptr<Particles3dcil<float>>& particles, std::vector<unsigned int> EmptyPlaces,
        int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<double>::InCell<Particles3dcil<double>>(
        const std::shared_ptr<Particles3dcil<double>>& particles, std::vector<unsigned int> EmptyPlaces,
        int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<float>::InCell<Particles2d<float>>(
        const std::shared_ptr<Particles2d<float>>& particles, std::vector<unsigned int> EmptyPlaces,
        int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<double>::InCell<Particles2d<double>>(
        const std::shared_ptr<Particles2d<double>>& particles, std::vector<unsigned int> EmptyPlaces,
        int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<float>::InCell<Particles2dpolar<float>>(
        const std::shared_ptr<Particles2dpolar<float>>& particles,
        std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<double>::InCell<Particles2dpolar<double>>(
        const std::shared_ptr<Particles2dpolar<double>>& particles,
        std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread, int flow);
template std::vector<unsigned int>
ParticleGridInterface<float>::InCellWithEps<Particles3dcil<float>>(
        const std::shared_ptr<Particles3dcil<float>>& particles, int i1, int i2, int thread, int flow);

template std::vector<unsigned int>
ParticleGridInterface<double>::InCellWithEps<Particles3dcil<double>>(
        const std::shared_ptr<Particles3dcil<double>>& particles, int i1, int i2, int thread, int flow);

template std::vector<unsigned int> ParticleGridInterface<float>::InCellWithEps<Particles2d<float>>(
        const std::shared_ptr<Particles2d<float>>& particles, int i1, int i2, int thread, int flow);

template std::vector<unsigned int>
ParticleGridInterface<double>::InCellWithEps<Particles2d<double>>(
        const std::shared_ptr<Particles2d<double>>& particles, int i1, int i2, int thread, int flow);

template std::vector<unsigned int>
ParticleGridInterface<float>::InCellWithEps<Particles2dpolar<float>>(
        const std::shared_ptr<Particles2dpolar<float>>& particles, int i1, int i2, int thread,
        int flow);

template std::vector<unsigned int>
ParticleGridInterface<double>::InCellWithEps<Particles2dpolar<double>>(
        const std::shared_ptr<Particles2dpolar<double>>& particles, int i1, int i2, int thread,
        int flow);

template void ParticleGridInterface<double>::axsPolar(
        Particles3dcil<double>* particles, const std::shared_ptr<GridData2dpolar<double>>& gridData,
        std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread);

template void ParticleGridInterface<double>::axsPolar(
        Particles3dcil<float>* particles, const std::shared_ptr<GridData2dpolar<double>>& gridData,
        std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread);