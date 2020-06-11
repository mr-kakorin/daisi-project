#include "EmissionCurrentSolver.h"
#include "Dmath.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice3d.h"
#include "Geom.h"
#include "GridData.h"
#include "ParticleGridInterface.h"
#include "ParticleSource.h"


template <class PointType>
template <class Archive>
void EmissionCurrentSolverBase<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& Lem;
    ar& Hem;
    ar& algorithm;
    ar& flowsNumbers;
}

template <class PointType>
template <class Archive>
void EmissionCurrentSolverBase<PointType>::load(Archive& ar, const unsigned int)
{
    ar& Lem;
    ar& Hem;
    ar& algorithm;
    ar& flowsNumbers;
}

template <class PointType>
void EmissionCurrentSolverBase<PointType>::CalculateCathodeFields(
    std::vector<std::shared_ptr<ParticleSource2d<PointType>>>& source,
    const std::shared_ptr<GridData3d<PointType>>& gridData, int flowNumber){

}

template <class PointType>
void EmissionCurrentSolverBase<PointType>::init(
    int i, const std::shared_ptr<GridData3d<PointType>>& gridData, int flagEm,
    std::vector<std::shared_ptr<ParticleSource2d<PointType>>> source,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface, int flag){

}

template <class PointType>
int EmissionCurrentSolverBase<PointType>::GetEmSize()
{
    int s = 0;
    for (int i = 0; i < nearCathodeVolumes.size(); i++)
    {
        if (0 == nearCathodeVolumes[i][0].flagCurrentLimited)
            s = s + nearCathodeVolumes[i].size();
    }
    return s;
}

template <class PointType>
void EmissionCurrentSolverBase<PointType>::SetParameters(std::vector<std::vector<double>> in)
{
    if (!in[0].size())
        return;

    algorithm = int(in[0][0]);
    for (int i = 0; i < Lem.size(); i++)
        Lem[i] = in[i + 1][0];
    for (int i = 0; i < Lem.size(); i++)
        Hem[i] = in[i + 1][1];
}

template <class PointType>
std::vector<std::vector<double>> EmissionCurrentSolverBase<PointType>::GetParameters()
{
    std::vector<std::vector<double>> result(4);
    result[0].push_back(int(algorithm));

    for (int i = 0; i < flowsNumbers.size(); i++)
        result[1].push_back(int(flowsNumbers[i]));

    result[2] = Lem;
    result[3] = Hem;

    return result;
}

template <class PointType>
void EmissionCurrentSolverBase<PointType>::addFlow(int flowNumber)
{
    Lem.push_back(0.001);
    Hem.push_back(0.001);
    flowsNumbers.push_back(flowNumber);
}

template <class PointType>
EmissionCurrentSolverBase<PointType>::EmissionCurrentSolverBase()
{
    algorithm = 0;
}

template <class PointType>
void EmissionCurrentSolverBase<PointType>::reset()
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
void EmissionCurrentSolverBase<PointType>::SetValueOnSource(
    const std::shared_ptr<ParticleSource2d<PointType>>& source, std::vector<double> value,
    int flowNumber, int flag)
{
    if (!source)
        return;
    if (value.size() == 1)
    {
        for (int kk = 0; kk < source->sourceSurface.size(); kk++)
        {
			//electric field
            if (flag == 0)
                source->sourceSurface[kk].E = value[0];
			//charge density
            if (flag == 1)
                source->sourceSurface[kk].currentDensity = value[0];
            //max charge density
			if (flag == 2)
                source->sourceSurface[kk].maximalCurrentDensity = value[0];
        }
        return;
    }

    for (int i = 0; i < points1[flowNumber].size() - 1; i++)
    {

        //	double Val = std::abs(value[i]);

        int dk = istarts[flowNumber][i + 1] - istarts[flowNumber][i];

        double Valold;
        double k0;
        double k1;
        int    kend;
        int    kstart;

        if (i >= 0)
        {
            k0 = (istarts[flowNumber][i] + istarts[flowNumber][i + 1]) / 2;
            k1 = (istarts[flowNumber][i + 1] + istarts[flowNumber][i + 2]) / 2;
            //	Valold = std::abs(value[i]);

            /*if (i == points1[flowNumber].size() - 1)
                    kend = istarts[flowNumber][i + 1];
            else
                    kend = k1;*/

            kend   = k1;
            kstart = k0;

            if (0 == i)
                kstart = 0;

            if (i == points1[flowNumber].size() - 2)
                kend = istarts[flowNumber][i + 2];

            // double Val1 = std::abs(value[i]);
            // double Val2 = std::abs(value[i+1]);

            double Val1 = value[i];
            double Val2 = value[i + 1];

            for (int kk = kstart; kk < kend; kk++)
            {
                //	double valInterp = (Val - Valold) * (kk - k0) / (k1 - k0) + Valold;
                double valInterp =
                    (Val1 * k1 - Val2 * k0) / (k1 - k0) + (Val1 - Val2) * kk / (k0 - k1);
                //	if (valInterp < 0)
                //	valInterp = 0.01;
                if (flag == 0)
                    source->sourceSurface[kk].E = valInterp;
                if (flag == 1)
                    source->sourceSurface[kk].currentDensity = valInterp;
                if (flag == 2)
                    source->sourceSurface[kk].maximalCurrentDensity = valInterp;
            }
        }
        /*if (i == 0)
        {
                for (int kk = 0; kk < istarts[flowNumber][i + 1] / 2; kk++)
                {
                        if (flag == 0)
                                source->sourceSurface[kk].E = Val;
                        if (flag == 1)
                                source->sourceSurface[kk].currentDensity = Val;
                }
        }*/
    }
}

template <class PointType>
template <class gridDataType>
void EmissionCurrentSolverBase<PointType>::CalculateCathodeFields(
    const std::shared_ptr<ParticleSource2d<PointType>>& source,
    const std::shared_ptr<gridDataType>& gridData, int flowNumber)
{
    double              Ex;
    double              Ey;
    double              ExCol;
    double              EyCol;
    double              ErAverage = 0;
    std::vector<double> test;
    for (int i = 0; i < points1[flowNumber].size(); i++)
    {
        gridData->interpolatePoint(nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
                                   nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Ex, Ey);
        //	gridData->interpolatePoint(nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
        // nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 1, ExCol, EyCol);

        CathodeFields[flowNumber][i] = (Ex * nearCathodeVolumes[flowNumber][i].normalX[0] +
                                        Ey * nearCathodeVolumes[flowNumber][i].normalY[0]);

        //	CathodeFields[flowNumber][i] = sqrt(Ex*Ex + Ey*Ey);
        ErAverage = ErAverage + std::abs(CathodeFields[flowNumber][i]);
    }
    ErAverage = ErAverage / points1[flowNumber].size();
    source->setErAverage(ErAverage);

    SetValueOnSource(source, CathodeFields[flowNumber], flowNumber, 0);
}

template <class PointType>
template <class gridDataType>
void EmissionCurrentSolverBase<PointType>::init(
    int i, const std::shared_ptr<gridDataType>& gridData, int flagEm,
    const std::shared_ptr<ParticleSource2d<PointType>>&      source,
    const std::shared_ptr<ParticleGridInterface<PointType>>& particleGridInterface, int flag)
{
    /*PointType r;
    PointType z;
    std::vector<unsigned int> indexesTmp(source->sourceSurface.size());
    for (int i = 0; i < source->sourceSurface.size(); i++)
    {
    r = source->sourceSurface[i].extractingEdge->Middle().x;
    z = source->sourceSurface[i].extractingEdge->Middle().y;
    indexesTmp[i] = particleGridInterface->InCell(r,z);
    }
    cellNumbers.push_back(indexesTmp);*/

    int flowNumber = -1;

    for (int j = 0; j < flowsNumbers.size(); j++)
    {
        if (flowsNumbers[j] == i)
            flowNumber = i;
    }

    if (flowNumber == -1)
        return;

    // int flowNumber = nearCathodeVolumes.size();

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

    double dL = Lem[flowNumber];

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
            edgeIntersect.point1 = source->sourceSurface[istartstmp.back()].extractingEdge->point1;
            edgeIntersect.point2 = source->sourceSurface[i].extractingEdge->point2;

            emittingCutEdge[flowNumber].push_back(edgeIntersect);

            int iaverage = ceil((i + istartstmp.back()) / 2);

            DGeo::Point<PointType> tmp = source->sourceSurface[iaverage].extractingEdge->Middle();

            //	if (flag == 3)
            //		Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);

            points1tmp.push_back(tmp);

            //		DGeo::Point<PointType> ptest = edgeIntersect.GetNormalPoint1(2 *
            // Hem[flowNumber]);

            //	if (flag == 3)
            //			Dmath::Cartesian2Polar(ptest.x, ptest.y, ptest.x, ptest.y);

            if (source->sourceSurface[iaverage].flagNormal == 1)
            {
                DGeo::Point<PointType> tmp = edgeIntersect.GetNormalPoint1(Hem[flowNumber]);
                nearCathodeVolumes[flowNumber].push_back(NearCathodeVolume<PointType>(
                    edgeIntersect, edgeIntersect.TranslateByNormal1(Hem[flowNumber]),
                    gridData->GetType(), flagEm));

                //		if (flag == 3)
                //			Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);
                points2tmp.push_back(tmp);
            }
            else
            {
                DGeo::Point<PointType> tmp = edgeIntersect.GetNormalPoint2(Hem[flowNumber]);

                nearCathodeVolumes[flowNumber].push_back(NearCathodeVolume<PointType>(
                    edgeIntersect, edgeIntersect.TranslateByNormal2(Hem[flowNumber]),
                    gridData->GetType(), flagEm));

                //	if (flag == 3)
                //		Dmath::Cartesian2Polar(tmp.x, tmp.y, tmp.x, tmp.y);
                points2tmp.push_back(tmp);
            }

            istartstmp.push_back(i + 1);
            n++;
        }
    }

    CathodeFields.back().resize(points1tmp.size());

    K.resize(flowNumber + 1);
    K.back().resize(points1tmp.size());
    for (int i      = 0; i < points1tmp.size(); i++)
        K.back()[i] = 0.1;

    istartstmp.back();
    gradients.emplace_back( DGeo::calc_grad2d( points1tmp ) );
    points1.push_back(points1tmp);
    points2.push_back(points2tmp);

    istarts.push_back(istartstmp);

    CalculateCathodeFields(source, gridData, flowNumber);
    E0.push_back(CathodeFields.back());
}

template class EmissionCurrentSolverBase<float>;
template class EmissionCurrentSolverBase<double>;

template void EmissionCurrentSolverBase<double>::CalculateCathodeFields<GridData2daxs<double>>(
        const std::shared_ptr<ParticleSource2d<double>>& source,
        const std::shared_ptr<GridData2daxs<double>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<float>::CalculateCathodeFields<GridData2daxs<float>>(
        const std::shared_ptr<ParticleSource2d<float>>& source,
        const std::shared_ptr<GridData2daxs<float>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<double>::CalculateCathodeFields<GridData2d<double>>(
        const std::shared_ptr<ParticleSource2d<double>>& source,
        const std::shared_ptr<GridData2d<double>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<float>::CalculateCathodeFields<GridData2d<float>>(
        const std::shared_ptr<ParticleSource2d<float>>& source,
        const std::shared_ptr<GridData2d<float>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<double>::CalculateCathodeFields<GridData2dpolar<double>>(
        const std::shared_ptr<ParticleSource2d<double>>& source,
        const std::shared_ptr<GridData2dpolar<double>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<float>::CalculateCathodeFields<GridData2dpolar<float>>(
        const std::shared_ptr<ParticleSource2d<float>>& source,
        const std::shared_ptr<GridData2dpolar<float>>& gridData, int flowNumber);
template void EmissionCurrentSolverBase<double>::init<GridData2daxs<double>>(
        int i, const std::shared_ptr<GridData2daxs<double>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<double>>&      source,
        const std::shared_ptr<ParticleGridInterface<double>>& particleGridInterface, int flag);
template void EmissionCurrentSolverBase<float>::init<GridData2daxs<float>>(
        int i, const std::shared_ptr<GridData2daxs<float>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<float>>&      source,
        const std::shared_ptr<ParticleGridInterface<float>>& particleGridInterface, int flag);
template void EmissionCurrentSolverBase<double>::init<GridData2d<double>>(
        int i, const std::shared_ptr<GridData2d<double>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<double>>&      source,
        const std::shared_ptr<ParticleGridInterface<double>>& particleGridInterface, int flag);
template void EmissionCurrentSolverBase<float>::init<GridData2d<float>>(
        int i, const std::shared_ptr<GridData2d<float>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<float>>&      source,
        const std::shared_ptr<ParticleGridInterface<float>>& particleGridInterface, int flag);
template void EmissionCurrentSolverBase<double>::init<GridData2dpolar<double>>(
        int i, const std::shared_ptr<GridData2dpolar<double>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<double>>&      source,
        const std::shared_ptr<ParticleGridInterface<double>>& particleGridInterface, int flag);
template void EmissionCurrentSolverBase<float>::init<GridData2dpolar<float>>(
        int i, const std::shared_ptr<GridData2dpolar<float>>& gridData, int flagEm,
        const std::shared_ptr<ParticleSource2d<float>>&      source,
        const std::shared_ptr<ParticleGridInterface<float>>& particleGridInterface, int flag);

template void EmissionCurrentSolverBase<float>::load<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void EmissionCurrentSolverBase<double>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void EmissionCurrentSolverBase<double>::load<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void EmissionCurrentSolverBase<float>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;