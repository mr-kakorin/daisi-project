#include "EmissionCurrentSolver.h"

template class EmissionCurrentSolver2dBase<float>;
template class EmissionCurrentSolver2dBase<double>;

template <class PointType>
void EmissionCurrentSolver2dBase<PointType>::reset()
{
    CathodeFields.clear();
    emissionCells.clear();
    cellNumbers.clear();
    points1.clear();
    points2.clear();
    istarts.clear();
    E0.clear();
    K.clear();
    nearCathodeVolumes.clear();
}

template <class PointType>
void EmissionCurrentSolver2dBase<PointType>::SetValueOnSource(ParticleSource2d<PointType>* source,
                                                              std::vector<double>          value,
                                                              int flowNumber, int flag)
{
    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        double Val = std::abs(value[i]);

        int dk = istarts[flowNumber][i + 1] - istarts[flowNumber][i];

        double Valold;
        double k0;
        double k1;

        if (i > 0)
        {
            k0     = (istarts[flowNumber][i] + istarts[flowNumber][i - 1]) / 2;
            k1     = (istarts[flowNumber][i + 1] + istarts[flowNumber][i]) / 2;
            Valold = std::abs(value[i - 1]);

            int kend;

            if (i == points1[flowNumber].size() - 1)
                kend = istarts[flowNumber][i + 1];
            else
                kend = k1;

            for (int kk = k0; kk < kend; kk++)
            {
                double valInterp = (Val - Valold) * (kk - k0) / (k1 - k0) + Valold;
                if (flag == 0)
                    source->sourceSurface[kk].E = valInterp;
                if (flag == 1)
                    source->sourceSurface[kk].currentDensity = valInterp;
            }
        }
        if (i == 0)
        {
            for (int kk = 0; kk < istarts[flowNumber][i + 1] / 2; kk++)
            {
                if (flag == 0)
                    source->sourceSurface[kk].E = Val;
                if (flag == 1)
                    source->sourceSurface[kk].currentDensity = Val;
            }
        }
    };
};

template void EmissionCurrentSolver2dBase<double>::CalculateCathodeFields<GridData2daxs<double>*>(
    ParticleSource2d<double>* source, GridData2daxs<double>* gridData, int flowNumber);
template void EmissionCurrentSolver2dBase<float>::CalculateCathodeFields<GridData2daxs<float>*>(
    ParticleSource2d<float>* source, GridData2daxs<float>* gridData, int flowNumber);
template void EmissionCurrentSolver2dBase<double>::CalculateCathodeFields<GridData2d<double>*>(
    ParticleSource2d<double>* source, GridData2d<double>* gridData, int flowNumber);
template void EmissionCurrentSolver2dBase<float>::CalculateCathodeFields<GridData2d<float>*>(
    ParticleSource2d<float>* source, GridData2d<float>* gridData, int flowNumber);
template void EmissionCurrentSolver2dBase<double>::CalculateCathodeFields<GridData2dpolar<double>*>(
    ParticleSource2d<double>* source, GridData2dpolar<double>* gridData, int flowNumber);
template void EmissionCurrentSolver2dBase<float>::CalculateCathodeFields<GridData2dpolar<float>*>(
    ParticleSource2d<float>* source, GridData2dpolar<float>* gridData, int flowNumber);

template <class PointType>
template <class gridDataType>
void EmissionCurrentSolver2dBase<PointType>::CalculateCathodeFields(
    ParticleSource2d<PointType>* source, gridDataType gridData, int flowNumber)
{
    double Ex;
    double Ey;
    double ErAverage = 0;
    for (int i = 0; i < points1[flowNumber].size(); i++)
    {
        gridData->interpolatePoint(nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
                                   nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Ex, Ey);
        CathodeFields[flowNumber][i] = (Ex * nearCathodeVolumes[flowNumber][i].normalX[0] +
                                        Ey * nearCathodeVolumes[flowNumber][i].normalY[0]);
        ErAverage                    = ErAverage + std::abs(CathodeFields[flowNumber][i]);
    }
    ErAverage = ErAverage / points1[flowNumber].size();
    source->setErAverage(ErAverage);

    SetValueOnSource(source, CathodeFields[flowNumber], flowNumber, 0);
};

template void EmissionCurrentSolver2dBase<double>::init<GridData2daxs<double>*>(
    GridData2daxs<double>* gridData, double dH, int flagEm, ParticleSource2d<double>* source,
    ParticleGridInterface2d<double>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);
template void EmissionCurrentSolver2dBase<float>::init<GridData2daxs<float>*>(
    GridData2daxs<float>* gridData, double dH, int flagEm, ParticleSource2d<float>* source,
    ParticleGridInterface2d<float>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);
template void EmissionCurrentSolver2dBase<double>::init<GridData2d<double>*>(
    GridData2d<double>* gridData, double dH, int flagEm, ParticleSource2d<double>* source,
    ParticleGridInterface2d<double>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);
template void EmissionCurrentSolver2dBase<float>::init<GridData2d<float>*>(
    GridData2d<float>* gridData, double dH, int flagEm, ParticleSource2d<float>* source,
    ParticleGridInterface2d<float>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);
template void EmissionCurrentSolver2dBase<double>::init<GridData2dpolar<double>*>(
    GridData2dpolar<double>* gridData, double dH, int flagEm, ParticleSource2d<double>* source,
    ParticleGridInterface2d<double>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);
template void EmissionCurrentSolver2dBase<float>::init<GridData2dpolar<float>*>(
    GridData2dpolar<float>* gridData, double dH, int flagEm, ParticleSource2d<float>* source,
    ParticleGridInterface2d<float>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells);

template <class PointType>
template <class gridDataType>
void EmissionCurrentSolver2dBase<PointType>::init(
    gridDataType gridData, double dH, int flagEm, ParticleSource2d<PointType>* source,
    ParticleGridInterface2d<PointType>* particleGridInterface, int flag,
    std::vector<unsigned int> realCells)
{
    /*PointType r;
    PointType z;
    std::vector<unsigned int> indexesTmp(source->sourceSurface.size());
    for (int i = 0; i < source->sourceSurface.size(); i++)
    {
    r = source->sourceSurface[i].extractingEdge.Middle().x;
    z = source->sourceSurface[i].extractingEdge.Middle().y;
    indexesTmp[i] = particleGridInterface->InCell(r,z);
    }
    cellNumbers.push_back(indexesTmp);*/

    int flowNumber = nearCathodeVolumes.size();

    nearCathodeVolumes.resize(flowNumber + 1);
    emittingCutEdge.resize(flowNumber + 1);
    CathodeFields.resize(flowNumber + 1);

    PointType r;
    PointType z;

    PointType dr;
    PointType dz;

    std::vector<DGeo::Point<PointType>> points1tmp;
    std::vector<DGeo::Point<PointType>> points2tmp;

    std::vector<int> istartstmp;

    std::vector<unsigned int> indexesTmp(source->sourceSurface.size());

    double L = source->sourceSurface.back().curveLength;

    double dL = L / 200.0;
    dL        = Lem;

    double LCurrent = 0;

    int n = 1;
    istartstmp.push_back(0);

    for (int i = 0; i < source->sourceSurface.size(); i++)
    {

        LCurrent = source->sourceSurface[i].curveLength;

        int flagLast = 1;
        if (istartstmp.size() > 1)
        {
            int lastdi = istartstmp[1] - istartstmp[0];
            if (source->sourceSurface.size() - i < 0.4 * lastdi)
                flagLast = 0;
        }

        if ((LCurrent >= dL * n && flagLast) || i == source->sourceSurface.size() - 1)
        {

            DGeo::Point<PointType> p1, p2;

            DGeo::Edge<PointType> edgeIntersect;
            edgeIntersect.point1 = source->sourceSurface[istartstmp.back()].extractingEdge.point1;
            edgeIntersect.point2 = source->sourceSurface[i].extractingEdge.point2;

            emittingCutEdge[flowNumber].push_back(edgeIntersect);

            int iaverage = ceil((i + istartstmp.back()) / 2);

            DGeo::Point<PointType> tmp = source->sourceSurface[iaverage].extractingEdge.Middle();

            if (flag == 3)
                Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);

            points1tmp.push_back(tmp);

            DGeo::Point<PointType> ptest = edgeIntersect.GetNormalPoint1(2 * dH);

            if (flag == 3)
                Dmath::Cartesian2Polar(ptest.x, ptest.y, ptest.x, ptest.y);

            if (particleGridInterface->InCellWithEps(ptest.x, ptest.y) != -1)
            {
                DGeo::Point<PointType> tmp = edgeIntersect.GetNormalPoint1(0.1 * dH);

                nearCathodeVolumes[flowNumber].push_back(NearCathodeVolume<PointType>(
                    edgeIntersect, edgeIntersect.TranslateByNormal1(Hem), gridData->GetType(),
                    flagEm));

                if (flag == 3)
                    Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);
                points2tmp.push_back(tmp);
            }
            else
            {
                DGeo::Point<PointType> tmp = edgeIntersect.GetNormalPoint2(0.1 * dH);

                nearCathodeVolumes[flowNumber].push_back(NearCathodeVolume<PointType>(
                    edgeIntersect, edgeIntersect.TranslateByNormal2(Hem), gridData->GetType(),
                    flagEm));

                if (flag == 3)
                    Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);
                points2tmp.push_back(tmp);
            }

            istartstmp.push_back(i + 1);
            n++;
        };
    }

    CathodeFields.back().resize(points1tmp.size());

    K.resize(flowNumber + 1);
    K.back().resize(points1tmp.size());
    for (int i = 0; i < points1tmp.size(); i++)
        K.back()[i] = 0.1;

    istartstmp.back();
    points1.push_back(points1tmp);
    points2.push_back(points2tmp);

    istarts.push_back(istartstmp);

    CalculateCathodeFields(source, gridData, flowNumber);
    E0.push_back(CathodeFields.back());
};
