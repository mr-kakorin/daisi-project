#include "Particle.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "GridData.h"
#include "ParticleGridInterface.h"
#include "ParticleShape2dCIC.h"
#include "ParticleShape2dTSC.h"
#include "WcalculateVector.h"


template <class PointType>
template <class particlesType>
std::vector<unsigned int>& ParticleGridInterface<PointType>::CheckParticlesBoundaries(
    const std::shared_ptr<BoundaryConditions>&                          boundaryConditions,
    std::vector<unsigned int>&                                          result,
    const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
    const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>&    conductorList,
    const std::shared_ptr<particlesType>& state1, const std::shared_ptr<particlesType>& state2,
    PointType dt, PointType charge, PointType mass, int i1, int i2, int thread)
{

    std::vector<int> list1 = boundaryConditions->GetDefaultConditionsList();
    int              s     = boundaryConditions->PropertyConditionListSize();

    std::vector<std::vector<int>> list;

    for (int i = 0; i < boundaryConditions->PropertyConditionListSize(); i++)
        list.push_back(boundaryConditions->GetPropertyConditionsBoundariesList(i));

    int                    tmp;
    DGeo::Point<PointType> tmpPoint;
    int                    flagContinue;

    std::vector<std::vector<unsigned int>>           particlesNumbers(s + 1);
    std::vector<std::vector<DGeo::Edge<PointType>>>  intersectionEdges(s + 1);
    std::vector<std::vector<DGeo::Point<PointType>>> intersectionPoints(s + 1);
    std::vector<std::vector<int>>                    intersectionBoundaryNumber(s + 1);

    std::vector<DGeo::Point<PointType>> intersectionPointsConductor;
    std::vector<unsigned int>           particlesNumbersConductor;
    std::vector<int>                    intersectionBoundaryNumberConductor;

    int                   flagNotBoundary;
    int                   indexLoc;
    bool                  flagIntersection;
    DGeo::Edge<PointType> traceEdge;

    std::vector<unsigned int> isEmitted;

    int boindaryChekLimit = 4;

    if (epsilon < 1e-2)
        boindaryChekLimit = 8;

    if (epsilon < 1e-5)
        boindaryChekLimit = 16;

    if (epsilon < 1e-8)
        boindaryChekLimit = 32;

    if (epsilon < 1e-11)
        boindaryChekLimit = 64;

    if (epsilon < 1e-14)
        boindaryChekLimit = 128;

    if (epsilon < 1e-17)
        boindaryChekLimit = 250;

    for (int k = i1; k < i2; k++)
    {
        if (state2->flagEmitted[k] < boindaryChekLimit + 2)
            state2->flagEmitted[k]++;

        if (int(state2->flagEmitted[k]) < boindaryChekLimit)
            continue;

        flagIntersection = false;
        indexLoc         = templNumb->getElem(state1->cellsNumbers[k - i1], 0, 0);

        if (indexLoc == -1)
            continue;

        DGeo::Point<PointType> ptm;
        flagContinue = 0;

        for (int i = 0; i < s; i++)
        {

            if (1 == boundaryConditions->GetPropertyConditionTypeFlag(i))
            {
                traceEdge.point1.x = state1->GetPointerToPosition1()[k - i1];
                traceEdge.point1.y = state1->GetPointerToPosition2()[k - i1];
                traceEdge.point1.z = 0;

                traceEdge.point2.x = state2->GetPointerToPosition1()[k];
                traceEdge.point2.y = state2->GetPointerToPosition2()[k];
                traceEdge.point2.z = 0;

                flagIntersection =
                    boundaryConditions->CheckBoundaryConditionIntersection(traceEdge, i);

                if (traceEdge.point2.x < 0.005)
                {
                    int tt = 0;
                }
                if (flagIntersection)
                {
                    particlesNumbers[i].push_back(k);
                    flagContinue = 1;

                    ptm.x = state2->GetPointerToCartesianX()[k];
                    ptm.y = state2->GetPointerToCartesianY()[k];

                    intersectionPoints[i].push_back(ptm);

                    intersectionBoundaryNumber[i].push_back(
                        boundaryConditions->GetElectrodeNumber(i));
                    break;
                }
            }
        }

        if (flagContinue)
            continue;

        switch (priorityParticleShapeType)
        {
        case 1:
            flagNotBoundary = TSCcellArray->isBoundary[indexLoc];
            break;
        case 0:
            flagNotBoundary = CICcellArray->isBoundary[indexLoc];
            break;
        }

        if (int(flagNotBoundary) == 0)
            continue;

        flagContinue       = 0;
        traceEdge.point1.x = state1->GetPointerToCartesianX()[k - i1];
        traceEdge.point1.y = state1->GetPointerToCartesianY()[k - i1];
        traceEdge.point1.z = 0;

        traceEdge.point2.x = state2->GetPointerToCartesianX()[k];
        traceEdge.point2.y = state2->GetPointerToCartesianY()[k];
        traceEdge.point2.z = 0;

        DGeo::Edge<PointType> EdgeTest;

        EdgeTest.point1.x = 0.029998030459810844;
        EdgeTest.point1.y = 0.070933231961481047;
        EdgeTest.point1.z = 0;

        EdgeTest.point2.x         = 0.030002170347779824;
        EdgeTest.point2.y         = 0.070935215495065637;
        EdgeTest.point2.z         = 0;
        bool flagIntersectionTest = false;

        for (int i = 0; i < s; i++)
        {
            for (int j = 0; j < list[i].size(); j++)
            {

                flagIntersection =
                    boundaries[list[i][j]]->IsIntersectionLight(epsilon, traceEdge, tmpPoint, tmp);

                if (flagIntersection)
                {
                    particlesNumbers[i].push_back(k);
                    intersectionEdges[i].push_back(boundaries[list[i][j]]->EdgesData[tmp]);
                    intersectionPoints[i].push_back(tmpPoint);
                    intersectionBoundaryNumber[i].push_back(-1);
                    flagContinue = 1;

                    for (int k = 0; k < conductorList.size(); k++)
                    {
                        for (int k1 = 0; k1 < conductorList[k]->boundaryNumbers.size(); k1++)
                        {
                            if (list[i][j] == conductorList[k]->boundaryNumbers[k1])
                            {
                                intersectionBoundaryNumber[i].back() = k;
                                break;
                            }
                        }
                    }
                    flagContinue = 1;
                    break;
                }
            }

            if (flagContinue)
                break;
        }

        if (flagContinue)
            continue;
    }

    for (int i = 0; i < s; i++)
    {
        ApplyBoundaryCondition(boundaryConditions, i, result,
                               boundaryConditions->GetConditionPropertyType(i), conductorList,
                               boundaryConditions->GetConditionPropertiesSimple(i), state1, state2,
                               particlesNumbers[i], intersectionEdges[i], intersectionPoints[i],
                               intersectionBoundaryNumber[i], dt, charge, mass, i1, thread);
    }

    for (int k = i1; k < i2; k++)
    {

        if (int(state2->flagEmitted[k]) < boindaryChekLimit)
        {
            continue;
        }

        flagContinue = 0;

        flagIntersection = false;
        indexLoc         = templNumb->getElem(state1->cellsNumbers[k - i1], 0, 0);

        if (indexLoc == -1)
            continue;

        switch (priorityParticleShapeType)
        {
        case 1:
            flagNotBoundary = TSCcellArray->isBoundary[indexLoc];
            break;
        case 0:
            flagNotBoundary = CICcellArray->isBoundary[indexLoc];
            break;
        }

        // if (int(flagNotBoundary) == 0 || state2->GetPointerToPosition1()[k]>0.101 ||
        // state2->GetPointerToPosition1()[k]<0.099)
        if (int(flagNotBoundary) == 0)
            continue;

        int flagC = 0;
        for (int tt = 0; tt < s; tt++)
        {
            for (int tt1 = 0; tt1 < particlesNumbers[tt].size(); tt1++)
            {
                if (particlesNumbers[tt][tt1] == k)
                {
                    flagC = 1;
                    break;
                }
            }
            if (flagC == 1)
                break;
        }

        if (flagC == 1)
            continue;

        traceEdge.point1.x = state1->GetPointerToCartesianX()[k - i1];
        traceEdge.point1.y = state1->GetPointerToCartesianY()[k - i1];
        traceEdge.point1.z = 0;

        traceEdge.point2.x = state2->GetPointerToCartesianX()[k];
        traceEdge.point2.y = state2->GetPointerToCartesianY()[k];
        traceEdge.point2.z = 0;

        if (k == 22)
        {
            int tt = 0;
        }

        for (int j = 0; j < list1.size(); j++)
        {

            flagIntersection =
                boundaries[list1[j]]->IsIntersectionLight(epsilon, traceEdge, tmpPoint, tmp);
            /*if (rt2 > rt1 && mass<1e-30)
            {
            flagIntersection = true;
            tmp = 0;
            tmpPoint = traceEdge.point2;
            }*/
            if (flagIntersection)
            {
                state2->SetCartesianPosition(k, {tmpPoint.x, tmpPoint.y});
                particlesNumbers[s].push_back(k);
                intersectionEdges[s].push_back(boundaries[list1[j]]->EdgesData[tmp]);
                break;
            }
        }

        if (flagIntersection == false)
            continue;
        flagContinue = 0;
        int flagTest = 0;
        for (int j = 0; j < conductorList.size(); j++)
        {
            for (int k1 = 0; k1 < conductorList[j]->boundaryNumbers.size(); k1++)
            {
                if (boundaries[conductorList[j]->boundaryNumbers[k1]]->IsIntersectionLight(
                        epsilon, traceEdge, tmpPoint, tmp))
                {
                    intersectionPointsConductor.push_back(tmpPoint);
                    particlesNumbersConductor.push_back(k);
                    intersectionBoundaryNumberConductor.push_back(j);
                    flagContinue = 1;
                    flagTest     = 1;
                    break;
                }
            }
            if (flagContinue)
                break;
        }
    }

    for (int j = 0; j < particlesNumbers[s].size(); j++)
    {
        int flag = 0;
        for (int i = 0; i < result.size(); i++)
        {
            if (result[i] == particlesNumbers[s][j])
                flag = 1;
        }
        if (flag == 0)
            result.push_back(particlesNumbers[s][j]);
    }

    //	ApplyDefaultCondition(conductorList, state1, state2, particlesNumbers[s],
    // intersectionEdges[s],
    // intersectionPoints[s], intersectionBoundaryNumber[s]);
    ApplyDefaultCondition(conductorList, state1, state2, particlesNumbersConductor,
                          intersectionEdges[s], intersectionPointsConductor,
                          intersectionBoundaryNumberConductor, dt, charge, mass, thread);

    /*for (int k = 0; k < state2->NParticles(); k++)
    {
            if (int(state2->flagEmitted[k]) == 2)
            {
                    state2->flagEmitted[k] = 1;
            }
    }*/

    return result;
}

template <class PointType>
void ParticleGridInterface<PointType>::SearchBoundariesCells(
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries)
{
    std::vector<unsigned int> searchIndexesLoc;

    searchIndexesLoc.resize(templNumb->ArSize());
    for (int k              = 0; k < templNumb->ArSize(); k++)
        searchIndexesLoc[k] = k;

    for (int j = 0; j < boundaries.size(); j++)
    {
        for (int i = 0; i < boundaries[j]->EdgesData.size(); i++)
        {
            DGeo::Point<PointType> p = boundaries[j]->EdgesData[i].Middle();

            int       cell = InCellWithEps(p.x, p.y, searchIndexesLoc);
            PointType H    = sqrt(GetH2(cell) * GetH2(cell) + GetH1(cell) * GetH1(cell));
        }
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::ApplyBoundaryCondition(
    const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
    std::vector<unsigned int>& empty, const std::string& conditionType,
    const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& conductorList,
    std::vector<double> conditionProperties, particlesType state1, particlesType state2,
    const std::vector<unsigned int>&           particlesNumbers,
    const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
    const std::vector<DGeo::Point<PointType>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
    PointType mass, int i1, int thread)
{
    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[0]) //"Metal Absorbtion"
    {
        if (state2->Get_currentCoef().size() == 0 || std::abs(dt - 1) > 0.000001)
        {
            for (int i = 0; i < particlesNumbers.size(); i++)
            {

                if (intersectionBoundaryNumber[i] != -1)
                    conductorList[intersectionBoundaryNumber[i]]->AddCharge(
                        state2->q[particlesNumbers[i]] * (1 - conditionProperties[0]),
                        intersectionPoints[i], dt, state2->GetEnergy(particlesNumbers[i], mass),
                        charge, thread);
                state2->q[particlesNumbers[i]] =
                    state2->q[particlesNumbers[i]] * conditionProperties[0];

                if (std::abs(state2->q[particlesNumbers[i]]) < fabsorpFrac * state2->avCharge ||
                    std::abs(state2->q[particlesNumbers[i]]) == 0)
                {
                    empty.push_back(particlesNumbers[i]);
                    state2->SetCartesianPosition(
                        particlesNumbers[i], {intersectionPoints[i].x, intersectionPoints[i].y});
                    state2->cellsNumbers[particlesNumbers[i]] = -1;
                }
            }
        }
        else
        {
            for (int i = 0; i < particlesNumbers.size(); i++)
            {

                if (intersectionBoundaryNumber[i] != -1)
                    conductorList[intersectionBoundaryNumber[i]]->AddCharge(
                        state2->q[particlesNumbers[i]] * (1 - conditionProperties[0]),
                        intersectionPoints[i], dt, state2->GetEnergy(particlesNumbers[i], mass),
                        charge, thread);
                state2->q[particlesNumbers[i]] =
                    state2->q[particlesNumbers[i]] * conditionProperties[0];
                state2->Get_currentCoef()[particlesNumbers[i]] =
                    state2->Get_currentCoef()[particlesNumbers[i]] * conditionProperties[0];

                if (std::abs(state2->q[particlesNumbers[i]]) < fabsorpFrac * state2->avCharge ||
                    std::abs(state2->q[particlesNumbers[i]]) == 0)
                {
                    empty.push_back(particlesNumbers[i]);
                    state2->SetCartesianPosition(
                        particlesNumbers[i], {intersectionPoints[i].x, intersectionPoints[i].y});
                    state2->cellsNumbers[particlesNumbers[i]] = -1;
                }
            }
        }
    }

    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[2]) //"backscattering"
    {
        for (int k = 0; k < particlesNumbers.size(); k++)
        {
            DGeo::Edge<PointType> traceEdge;
            traceEdge.point1.x = state1->GetPointerToCartesianX()[particlesNumbers[k] - i1];
            traceEdge.point1.y = state1->GetPointerToCartesianY()[particlesNumbers[k] - i1];
            traceEdge.point1.z = 0;

            traceEdge.point2.x   = state2->GetPointerToCartesianX()[particlesNumbers[k]];
            traceEdge.point2.y   = state2->GetPointerToCartesianY()[particlesNumbers[k]];
            traceEdge.point2.z   = 0;
            PointType alphaTrace = traceEdge.alpha();
            PointType alphaEdge  = intersectionEdges[k].alpha();
            PointType dalpha     = alphaEdge - alphaTrace;

            DGeo::Point<PointType> newPoint2 =
                traceEdge.point2.rotate(2 * dalpha, intersectionPoints[k]);
            state2->SetCartesianPosition(particlesNumbers[k], {newPoint2.x, newPoint2.y});

            DGeo::Point<PointType> velocityPoint1 = traceEdge.Middle();
            DGeo::Point<PointType> velocityPoint2;
            velocityPoint2.x = velocityPoint1.x + state2->GetCartesianPX(particlesNumbers[k]);
            velocityPoint2.y = velocityPoint1.y + state2->GetCartesianPY(particlesNumbers[k]);
            velocityPoint2.z = 0;

            DGeo::Point<PointType> newvelocityPoint1 =
                velocityPoint1.rotate(2 * dalpha, intersectionPoints[k]);
            DGeo::Point<PointType> newvelocityPoint2 =
                velocityPoint2.rotate(2 * dalpha, intersectionPoints[k]);

            state2->SetCartesianMomentum(particlesNumbers[k],
                                         {newvelocityPoint2.x - newvelocityPoint1.x,
                                          newvelocityPoint2.y - newvelocityPoint1.y});

            if (intersectionBoundaryNumber[k] != -1)
                conductorList[intersectionBoundaryNumber[k]]->AddCharge(
                    state2->q[particlesNumbers[k]] * (1 - conditionProperties[0]),
                    intersectionPoints[k], dt, state2->GetEnergy(particlesNumbers[k], mass), charge,
                    thread);

            state2->MultiplyMomentum(particlesNumbers[k], sqrt(conditionProperties[1]));

            state2->q[particlesNumbers[k]] =
                state2->q[particlesNumbers[k]] * conditionProperties[0];

            if (!(state2->Get_currentCoef().size() == 0 || std::abs(dt - 1) > 0.000001))
                state2->Get_currentCoef()[particlesNumbers[k]] =
                    state2->Get_currentCoef()[particlesNumbers[k]] * conditionProperties[0];

            if (std::abs(state2->q[particlesNumbers[k]]) < fabsorpFrac * state2->avCharge ||
                std::abs(state2->q[particlesNumbers[k]]) == 0)
            {
                empty.push_back(particlesNumbers[k]);
                state2->SetCartesianPosition(particlesNumbers[k],
                                             {intersectionPoints[k].x, intersectionPoints[k].y});
            }
        }
    }

    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[1]) // "Reflection"
    {
        if (1 == boundaryConditions->GetPropertyConditionTypeFlag(flag))
        {
            for (int k = 0; k < particlesNumbers.size(); k++)
            {
                state2->GetPointerToPosition1()[particlesNumbers[k]] =
                    state1->GetPointerToPosition1()[particlesNumbers[k] - i1];
                state2->GetPointerToPosition2()[particlesNumbers[k]] =
                    state1->GetPointerToPosition2()[particlesNumbers[k] - i1];
            }
            if (0 == int(boundaryConditions->GetPropertyConditionManualRestictions(flag)[0]))
            {
                for (int k = 0; k < particlesNumbers.size(); k++)
                {
                    state2->GetPointerToMomentum1()[particlesNumbers[k]] =
                        -state2->GetPointerToMomentum1()[particlesNumbers[k]];
                    state2->flagEmitted[particlesNumbers[k]] = 0;
                }
            }
            if (1 == int(boundaryConditions->GetPropertyConditionManualRestictions(flag)[0]))
            {
                for (int k = 0; k < particlesNumbers.size(); k++)
                {
                    state2->GetPointerToMomentum2()[particlesNumbers[k]] =
                        -state2->GetPointerToMomentum2()[particlesNumbers[k]];
                    state2->flagEmitted[particlesNumbers[k]] = 0;
                }
            }
        }
        else
        {
            for (int k = 0; k < particlesNumbers.size(); k++)
            {
                DGeo::Edge<PointType> traceEdge;
                traceEdge.point1.x = state1->GetPointerToCartesianX()[particlesNumbers[k] - i1];
                traceEdge.point1.y = state1->GetPointerToCartesianY()[particlesNumbers[k] - i1];
                traceEdge.point1.z = 0;

                traceEdge.point2.x   = state2->GetPointerToCartesianX()[particlesNumbers[k]];
                traceEdge.point2.y   = state2->GetPointerToCartesianY()[particlesNumbers[k]];
                traceEdge.point2.z   = 0;
                PointType alphaTrace = traceEdge.alpha();
                PointType alphaEdge  = intersectionEdges[k].alpha();
                PointType dalpha     = alphaEdge - alphaTrace;

                DGeo::Point<PointType> newPoint2 =
                    traceEdge.point2.rotate(2 * dalpha, intersectionPoints[k]);
                //	newPoint2 = traceEdge.point1;
                state2->SetCartesianPosition(particlesNumbers[k], {newPoint2.x, newPoint2.y});

                DGeo::Point<PointType> velocityPoint1 = traceEdge.Middle();
                DGeo::Point<PointType> velocityPoint2;

                PointType px = state2->GetCartesianPX(particlesNumbers[k]);
                PointType py = state2->GetCartesianPY(particlesNumbers[k]);

                velocityPoint2.x = velocityPoint1.x + state2->GetCartesianPX(particlesNumbers[k]);
                velocityPoint2.y = velocityPoint1.y + state2->GetCartesianPY(particlesNumbers[k]);
                velocityPoint2.z = 0;

                DGeo::Point<PointType> newvelocityPoint1 =
                    velocityPoint1.rotate(2 * dalpha, intersectionPoints[k]);
                DGeo::Point<PointType> newvelocityPoint2 =
                    velocityPoint2.rotate(2 * dalpha, intersectionPoints[k]);

                state2->SetCartesianMomentum(particlesNumbers[k],
                                             {newvelocityPoint2.x - newvelocityPoint1.x,
                                              newvelocityPoint2.y - newvelocityPoint1.y});
                //	state2->SetCartesianMomentum(particlesNumbers[k],
                // state2->GetCartesianPX(particlesNumbers[k],
                // newvelocityPoint2.y - newvelocityPoint1.y));
                // state2->GetPointerToMomentum1()[particlesNumbers[k]] =
                //-state2->GetPointerToMomentum1()[particlesNumbers[k]];
                state2->flagEmitted[particlesNumbers[k]] = 0;
            }
        }
    }
}

template <class PointType>
template <class particlesType>
void ParticleGridInterface<PointType>::ApplyDefaultCondition(
    const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& conductorList,
    particlesType state1, particlesType state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
    const std::vector<DGeo::Point<PointType>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
    PointType mass, int thread)
{
    for (int i = 0; i < particlesNumbers.size(); i++)
    {
        state2->cellsNumbers[particlesNumbers[i]] = -1;
        if (intersectionBoundaryNumber[i] != -1)
            conductorList[intersectionBoundaryNumber[i]]->AddCharge(
                state2->q[particlesNumbers[i]], intersectionPoints[i], dt,
                state2->GetEnergy(particlesNumbers[i], mass), charge, thread);
    }
    //	state2->EmptyPlaces = particlesNumbers;
}

template class ParticleGridInterface<float>;
template class ParticleGridInterface<double>;

template std::vector<unsigned int>&
ParticleGridInterface<float>::CheckParticlesBoundaries<Particles3dcil<float>>(
        const std::shared_ptr<BoundaryConditions>&                      boundaryConditions,
        std::vector<unsigned int>&                                      result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<float>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>&    conductorList,
        const std::shared_ptr<Particles3dcil<float>>&                   state1,
        const std::shared_ptr<Particles3dcil<float>>& state2, float dt, float charge, float mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<float>::CheckParticlesBoundaries<Particles2d<float>>(
        const std::shared_ptr<BoundaryConditions>&                      boundaryConditions,
        std::vector<unsigned int>&                                      result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<float>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>&    conductorList,
        const std::shared_ptr<Particles2d<float>>&                      state1,
        const std::shared_ptr<Particles2d<float>>& state2, float dt, float charge, float mass, int i1,
        int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<float>::CheckParticlesBoundaries<Particles2dpolar<float>>(
        const std::shared_ptr<BoundaryConditions>&                      boundaryConditions,
        std::vector<unsigned int>&                                      result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<float>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>&    conductorList,
        const std::shared_ptr<Particles2dpolar<float>>&                 state1,
        const std::shared_ptr<Particles2dpolar<float>>& state2, float dt, float charge, float mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<double>::CheckParticlesBoundaries<Particles3dcil<double>>(
        const std::shared_ptr<BoundaryConditions>&                       boundaryConditions,
        std::vector<unsigned int>&                                       result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<double>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>&    conductorList,
        const std::shared_ptr<Particles3dcil<double>>&                   state1,
        const std::shared_ptr<Particles3dcil<double>>& state2, double dt, double charge, double mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<double>::CheckParticlesBoundaries<Particles2d<double>>(
        const std::shared_ptr<BoundaryConditions>&                       boundaryConditions,
        std::vector<unsigned int>&                                       result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<double>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>&    conductorList,
        const std::shared_ptr<Particles2d<double>>&                      state1,
        const std::shared_ptr<Particles2d<double>>& state2, double dt, double charge, double mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<double>::CheckParticlesBoundaries<Particles2dpolar<double>>(
        const std::shared_ptr<BoundaryConditions>&                       boundaryConditions,
        std::vector<unsigned int>&                                       result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<double>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>&    conductorList,
        const std::shared_ptr<Particles2dpolar<double>>&                 state1,
        const std::shared_ptr<Particles2dpolar<double>>& state2, double dt, double charge, double mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<double>::CheckParticlesBoundaries<Particles3d<double>>(
        const std::shared_ptr<BoundaryConditions>&                       boundaryConditions,
        std::vector<unsigned int>&                                       result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<double>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>&    conductorList,
        const std::shared_ptr<Particles3d<double>>&                      state1,
        const std::shared_ptr<Particles3d<double>>& state2, double dt, double charge, double mass,
        int i1, int i2, int thread);

template std::vector<unsigned int>&
ParticleGridInterface<float>::CheckParticlesBoundaries<Particles3d<float>>(
        const std::shared_ptr<BoundaryConditions>&                      boundaryConditions,
        std::vector<unsigned int>&                                      result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<float>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>&    conductorList,
        const std::shared_ptr<Particles3d<float>>&                      state1,
        const std::shared_ptr<Particles3d<float>>& state2, float dt, float charge, float mass, int i1,
        int i2, int thread);

template void ParticleGridInterface<float>::ApplyDefaultCondition<Particles3dcil<float>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        Particles3dcil<float>* state1, Particles3dcil<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass,
        int thread);

template void ParticleGridInterface<float>::ApplyDefaultCondition<Particles2d<float>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        Particles2d<float>* state1, Particles2d<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass,
        int thread);

template void ParticleGridInterface<float>::ApplyDefaultCondition<Particles2dpolar<float>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        Particles2dpolar<float>* state1, Particles2dpolar<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass,
        int thread);

template void ParticleGridInterface<double>::ApplyDefaultCondition<Particles3dcil<double>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        Particles3dcil<double>* state1, Particles3dcil<double>* state2,
        const std::vector<unsigned int>&        particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int thread);

template void ParticleGridInterface<double>::ApplyDefaultCondition<Particles2d<double>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        Particles2d<double>* state1, Particles2d<double>* state2,
        const std::vector<unsigned int>&        particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int thread);

template void ParticleGridInterface<double>::ApplyDefaultCondition<Particles2dpolar<double>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        Particles2dpolar<double>* state1, Particles2dpolar<double>* state2,
        const std::vector<unsigned int>&        particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int thread);

template void ParticleGridInterface<double>::ApplyDefaultCondition<Particles3d<double>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        Particles3d<double>* state1, Particles3d<double>* state2,
        const std::vector<unsigned int>&        particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int thread);

template void ParticleGridInterface<float>::ApplyDefaultCondition<Particles3d<float>*>(
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        Particles3d<float>* state1, Particles3d<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass,
        int thread);


template void ParticleGridInterface<float>::ApplyBoundaryCondition<Particles3dcil<float>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        std::vector<double> conditionProperties, Particles3dcil<float>* state1,
        Particles3dcil<float>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass, int i1,
        int thread);

template void ParticleGridInterface<float>::ApplyBoundaryCondition<Particles2d<float>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        std::vector<double> conditionProperties, Particles2d<float>* state1, Particles2d<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass, int i1,
        int thread);

template void ParticleGridInterface<float>::ApplyBoundaryCondition<Particles2dpolar<float>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        std::vector<double> conditionProperties, Particles2dpolar<float>* state1,
        Particles2dpolar<float>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass, int i1,
        int thread);

template void ParticleGridInterface<double>::ApplyBoundaryCondition<Particles3dcil<double>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        std::vector<double> conditionProperties, Particles3dcil<double>* state1,
        Particles3dcil<double>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int i1, int thread);

template void ParticleGridInterface<double>::ApplyBoundaryCondition<Particles2d<double>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        std::vector<double> conditionProperties, Particles2d<double>* state1,
        Particles2d<double>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int i1, int thread);

template void ParticleGridInterface<double>::ApplyBoundaryCondition<Particles2dpolar<double>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        std::vector<double> conditionProperties, Particles2dpolar<double>* state1,
        Particles2dpolar<double>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int i1, int thread);

template void ParticleGridInterface<double>::ApplyBoundaryCondition<Particles3d<double>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<double>>>& conductorList,
        std::vector<double> conditionProperties, Particles3d<double>* state1,
        Particles3d<double>* state2, const std::vector<unsigned int>& particlesNumbers,
        const std::vector<DGeo::Edge<double>>&  intersectionEdges,
        const std::vector<DGeo::Point<double>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, double dt, double charge, double mass,
        int i1, int thread);

template void ParticleGridInterface<float>::ApplyBoundaryCondition<Particles3d<float>*>(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<float>>>& conductorList,
        std::vector<double> conditionProperties, Particles3d<float>* state1, Particles3d<float>* state2,
        const std::vector<unsigned int>&       particlesNumbers,
        const std::vector<DGeo::Edge<float>>&  intersectionEdges,
        const std::vector<DGeo::Point<float>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass, int i1,
        int thread);