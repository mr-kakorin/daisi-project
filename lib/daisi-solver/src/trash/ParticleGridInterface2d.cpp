#include "ParticleGridInterface2d.h"
#include "WcalculateVector.h"

template ParticleGridInterface2d<float>;
template ParticleGridInterface2d<double>;

bool isTemplateEqual(std::vector<int> templ1, std::vector<int> templ2)
{
    for (int i = 0; i < 8; i++)
    {
        if ((templ1[i] >= 0 && templ2[i] < 0) || (templ1[i] < 0 && templ2[i] >= 0))
            return false;
    }
    return true;
};

template void ParticleGridInterface2d<float>::init<GridData2d<float>*>(GridData2d<float>* gridData,
                                                                       Dmath::imat templNumbIn, Dmath::imat flagMatrix,
                                                                       const std::vector<int>&           boundarypoints,
                                                                       const BoundaryContainer2d<float>& domainBoundary,
                                                                       int                               size);

template void ParticleGridInterface2d<float>::init<GridData2daxs<float>*>(
    GridData2daxs<float>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<float>& domainBoundary, int size);

template void ParticleGridInterface2d<float>::init<GridData2dpolar<float>*>(
    GridData2dpolar<float>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<float>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2d<double>*>(
    GridData2d<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2daxs<double>*>(
    GridData2daxs<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2dpolar<double>*>(
    GridData2dpolar<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template void ParticleGridInterface2d<float>::init<GridData2d<float>*>(GridData2d<float>* gridData,
                                                                       Dmath::imat templNumbIn, Dmath::imat flagMatrix,
                                                                       const std::vector<int>&           boundarypoints,
                                                                       const BoundaryContainer2d<float>& domainBoundary,
                                                                       int                               size);

template void ParticleGridInterface2d<float>::init<GridData2daxs<float>*>(
    GridData2daxs<float>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<float>& domainBoundary, int size);

template void ParticleGridInterface2d<float>::init<GridData2dpolar<float>*>(
    GridData2dpolar<float>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<float>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2d<double>*>(
    GridData2d<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2daxs<double>*>(
    GridData2daxs<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template void ParticleGridInterface2d<double>::init<GridData2dpolar<double>*>(
    GridData2dpolar<double>* gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
    const std::vector<int>& boundarypoints, const BoundaryContainer2d<double>& domainBoundary, int size);

template <class PointType>
template <class gridDatatype>
void ParticleGridInterface2d<PointType>::init(gridDatatype gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
                                              const std::vector<int>&               boundarypoints,
                                              const BoundaryContainer2d<PointType>& domainBoundary, int size)
{

    indexes.reserve(size);
    cells.reserve(size);

    searchIndexes.resize(9);

    x1Array = gridData->GetX1Pointer();
    x2Array = gridData->GetX2Pointer();

    templNumb = templNumbIn;
    std::vector<std::vector<int>> templCICBoundary;
    std::vector<int>              vtmp;

    vtmp = {-1, -1, -1, -1, 0, 0, -1, 0, 0};
    templCICBoundary.push_back(vtmp);
    vtmp = {-1, -1, -1, 0, 0, -1, 0, 0, -1};
    templCICBoundary.push_back(vtmp);
    vtmp = {0, 0, -1, 0, 0, -1, -1, -1, -1};
    templCICBoundary.push_back(vtmp);
    vtmp = {-1, 0, 0, -1, 0, 0, -1, -1, -1};
    templCICBoundary.push_back(vtmp);

    std::vector<int> templTSC = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<std::vector<int>> templCICx2TSCx1Boundary;
    vtmp = {0, 0, 0, 0, 0, 0, -1, -1, -1};
    templCICx2TSCx1Boundary.push_back(vtmp);
    vtmp = {-1, -1, -1, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, -1, -1, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {-1, -1, 0, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, 0, 0, -1, -1};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, 0, -1, -1, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    std::vector<std::vector<int>> templCICx1TSCx2Boundary;
    vtmp = {0, 0, -1, 0, 0, -1, 0, 0, -1};
    templCICx1TSCx2Boundary.push_back(vtmp);
    vtmp = {-1, 0, 0, -1, 0, 0, -1, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, -1, 0, 0, -1, 0, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, -1, 0, 0, -1};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {-1, 0, 0, -1, 0, 0, 0, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, -1, 0, 0, -1, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    std::vector<std::vector<int>> templ3CellBoundary;
    vtmp = {0, 0, 0, 0, 0, 0, 0, 0, -1};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {0, 0, 0, 0, 0, 0, -1, 0, 0};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {-1, 0, 0, 0, 0, 0, 0, 0, 0};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {0, 0, -1, 0, 0, 0, 0, 0, 0};
    templ3CellBoundary.push_back(vtmp);

    int n         = 0;
    int levelLow  = 0;
    int levelHigh = 0;

    cellArray.clear();
    iter = 0;
    CICArray.clear();

    int flagBreak;

    std::vector<int> flagOut;
    for (int j = 1; j < templNumb.GetNcol() - 1; j++)
    {

        for (int i = 1; i < templNumb.GetNrow() - 1; i++)
        {
            if (templNumb(i, j) != -1)
                flagOut.push_back(flagMatrix(i, j));
        };
    };

    switch (priorityParticleShapeType)
    {
    case 0:
        CICcellArray.init(flagOut.size(), gridData->GetX1Pointer(), gridData->GetX2Pointer());
        for (int j = 1; j < templNumb.GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb.GetNrow() - 1; i++)
            {
                if (templNumb(i, j) != -1)
                {
                    std::vector<int> v(9);
                    v[0] = templNumb(i - 1, j - 1);
                    v[1] = templNumb(i, j - 1);
                    v[2] = templNumb(i + 1, j - 1);
                    v[3] = templNumb(i - 1, j);
                    v[4] = templNumb(i, j);
                    v[5] = templNumb(i + 1, j);
                    v[6] = templNumb(i - 1, j + 1);
                    v[7] = templNumb(i, j + 1);
                    v[8] = templNumb(i + 1, j + 1);

                    levelHigh = templNumb(i, j + 1);
                    CICcellArray.AddCell(gridData, levelHigh, flagOut, v, domainBoundary, boundarypoints, 1);
                };
                n++;
            };
        };
        break;
    case 1:
        TSCcellArray.init(flagOut.size(), gridData->GetX1Pointer(), gridData->GetX2Pointer());
        for (int j = 1; j < templNumb.GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb.GetNrow() - 1; i++)
            {

                if (templNumb(i, j) != -1)
                {
                    flagBreak = 0;

                    std::vector<int> v(9);
                    v[0]                  = templNumb(i - 1, j - 1);
                    v[1]                  = templNumb(i, j - 1);
                    v[2]                  = templNumb(i + 1, j - 1);
                    v[3]                  = templNumb(i - 1, j);
                    v[4]                  = templNumb(i, j);
                    v[5]                  = templNumb(i + 1, j);
                    v[6]                  = templNumb(i - 1, j + 1);
                    v[7]                  = templNumb(i, j + 1);
                    v[8]                  = templNumb(i + 1, j + 1);
                    std::vector<int> vtmp = v;
                    int              flag = 1;

                    levelHigh = templNumb(i, j + 1);
                    levelLow  = templNumb(i, j - 1);

                    while (1)
                    {

                        if (isTemplateEqual(v, templTSC))
                        {
                            std::vector<int> t;
                            TSCcellArray.AddCell(gridData, levelLow, levelHigh, 0, t, flagOut, domainBoundary,
                                                 boundarypoints, 1);
                            // AddCell(gridData, 1, templNumb(i, j), levelLow, levelHigh);
                            flag = 0;
                            break;
                        };

                        if (isTemplateEqual(v, templCICBoundary[0]) || isTemplateEqual(v, templCICBoundary[1]) ||
                            isTemplateEqual(v, templCICBoundary[2]) || isTemplateEqual(v, templCICBoundary[3]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray.AddCell(gridData, levelLow, levelHigh, 1, v, flagOut, domainBoundary,
                                                 boundarypoints, 1);
                            // AddCell(gridData, templNumb(i, j), v[0], v[1], v[2], v[3]);
                            flag = 0;
                            break;
                        };

                        if (isTemplateEqual(v, templCICx2TSCx1Boundary[0]) ||
                            isTemplateEqual(v, templCICx2TSCx1Boundary[1]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray.AddCell(gridData, levelLow, levelHigh, 2, v, flagOut, domainBoundary,
                                                 boundarypoints, 1);
                            // AddCell(gridData, 1, templNumb(i, j), v[0], v[1], v[2], v[3], v[4], v[5]);
                            flag = 0;
                            break;
                        };

                        for (int k = 2; k < 6; k++)
                        {
                            if (isTemplateEqual(v, templCICx2TSCx1Boundary[k]))
                            {
                                switch (k)
                                {
                                case 2:
                                    v.erase(v.begin() + 0);
                                    break;
                                case 3:
                                    v.erase(v.begin() + 2);
                                    break;
                                case 4:
                                    v.erase(v.begin() + 6);
                                    break;
                                case 5:
                                    v.erase(v.begin() + 8);
                                    break;
                                }
                                size_t n = v.size();
                                for (int k = 0; k < n; k++)
                                {
                                    if (v[k] < 0)
                                    {
                                        v.erase(v.begin() + k);
                                        k--;
                                        n = v.size();
                                    };
                                }
                                TSCcellArray.AddCell(gridData, levelLow, levelHigh, 2, v, flagOut, domainBoundary,
                                                     boundarypoints, 1);
                                //	AddCell(gridData, 1, templNumb(i, j), v[0], v[1], v[2], v[3], v[4], v[5]);
                                flag      = 0;
                                flagBreak = 1;
                                break;
                            };
                        }
                        if (flagBreak)
                            break;

                        if (isTemplateEqual(v, templCICx1TSCx2Boundary[0]) ||
                            isTemplateEqual(v, templCICx1TSCx2Boundary[1]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray.AddCell(gridData, levelLow, levelHigh, 3, v, flagOut, domainBoundary,
                                                 boundarypoints, 1);
                            // AddCell(gridData, 0, templNumb(i, j), v[0], v[1], v[2], v[3], v[4], v[5]);
                            flag = 0;
                            break;
                        };

                        for (int k = 2; k < 6; k++)
                        {
                            if (isTemplateEqual(v, templCICx1TSCx2Boundary[k]))
                            {
                                switch (k)
                                {
                                case 2:
                                    v.erase(v.begin() + 8);
                                    break;
                                case 3:
                                    v.erase(v.begin() + 2);
                                    break;
                                case 4:
                                    v.erase(v.begin() + 6);
                                    break;
                                case 5:
                                    v.erase(v.begin() + 0);
                                    break;
                                }
                                size_t n = v.size();
                                for (int k = 0; k < n; k++)
                                {
                                    if (v[k] < 0)
                                    {
                                        v.erase(v.begin() + k);
                                        k--;
                                        n = v.size();
                                    };
                                }
                                TSCcellArray.AddCell(gridData, levelLow, levelHigh, 3, v, flagOut, domainBoundary,
                                                     boundarypoints, 1);
                                // AddCell(gridData, 0, templNumb(i, j), v[0], v[1], v[2], v[3], v[4], v[5]);
                                flag      = 0;
                                flagBreak = 1;
                                break;
                            };
                        }
                        if (flagBreak)
                            break;

                        if (isTemplateEqual(v, templ3CellBoundary[0]) || isTemplateEqual(v, templ3CellBoundary[1]) ||
                            isTemplateEqual(v, templ3CellBoundary[2]) || isTemplateEqual(v, templ3CellBoundary[3]))
                        {
                            TSCcellArray.AddCell(gridData, levelLow, levelHigh, 4, v, flagOut, domainBoundary,
                                                 boundarypoints, 1);
                            // AddCell(gridData, templNumb(i, j), v);
                            flag = 0;
                            break;
                        };
                    }
                    if (flag == 0)
                    {
                        int flag1 = 0;
                        for (int s = 0; s < boundarypoints.size(); s++)
                        {
                            for (int s1 = 0; s1 < 9; s1++)
                            {
                                if (boundarypoints[s] == vtmp[s1])
                                    flag1 = 1;
                            }
                        }
                        if (flag1 == 1)
                            isBoundary.push_back(true);
                        else
                            isBoundary.push_back(false);
                    };
                };
                n++;
            };
        };
        break;
    }
    // particleGridInterface->AddCell(0, 0, 0, 0);
    // particleGridInterface->AddCell(1, 0, 0, 0);
    //	ParticleGridInterface->cellArray[0]->InCell(0, 0);
    //	ParticleGridInterface->cellArray[1]->InCell(0, 0);
};

template <class PointType> int ParticleGridInterface2d<PointType>::InCell(PointType x1, PointType x2)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray.InCell(x1, x2);
    case 0:
        return CICcellArray.InCell(x1, x2);
    }
    return -1;
};

template <class PointType> int ParticleGridInterface2d<PointType>::InCellWithEps(PointType x1, PointType x2)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray.InCellWithEps(x1, x2);
    case 0:
        return CICcellArray.InCellWithEps(x1, x2);
    }
    return -1;
};

template <class PointType> PointType ParticleGridInterface2d<PointType>::GetH1(int index)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray.GetH1(index);
    case 0:
        return CICcellArray.GetH1(index);
    }
};

template <class PointType> PointType ParticleGridInterface2d<PointType>::GetH2(int index)
{
    switch (priorityParticleShapeType)
    {
    case 1:
        return TSCcellArray.GetH2(index);
    case 0:
        return CICcellArray.GetH2(index);
    }
};

template <class PointType>
int ParticleGridInterface2d<PointType>::InCell(PointType x1, PointType x2, const std::vector<int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb.getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (TSCcellArray.InCell(index, x1, x2))
                    return searchIndexes[i];
            }
        };
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb.getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (CICcellArray.InCell(index, x1, x2))
                    return searchIndexes[i];
            }
        };
        break;
    }
    return -1;
};

template <class PointType>
int ParticleGridInterface2d<PointType>::InCellWithEps(PointType x1, PointType x2, const std::vector<int>& searchIndexes)
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb.getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (TSCcellArray.InCellWithEps(index, x1, x2))
                    return searchIndexes[i];
            }
        };
        break;
    case 0:
        for (int i = 0; i < searchIndexes.size(); i++)
        {
            index = templNumb.getElem(searchIndexes[i], 0, 0);
            if (index != -1)
            {
                if (CICcellArray.InCellWithEps(index, x1, x2))
                    return searchIndexes[i];
            }
        };
        break;
    }
    return -1;
};

/*int ParticleGridInterface2d <PointType>::InCellNew(PointType x1, PointType x2, const std::vector<int>& searchIndexes,
const GridData2daxs<DataContainer>* const gridData)
{
        switch (priorityParticleShapeType)
        {
        case 1:
                for (int i = 0; i < searchIndexes.size(); i++)
                {
                        int index = templNumb.getElem(searchIndexes[i], 0, 0);
                        if (index != -1)
                        {
                                if (cellArray[index]->InCell(x1, x2, gridData->r, gridData->z))
                                        return searchIndexes[i];
                        }
                };
                break;
        case 0:
                for (int i = 0; i < searchIndexes.size(); i++)
                {
                        int index = templNumb.getElem(searchIndexes[i], 0, 0);
                        if (index != -1)
                        {
                                if (ParticleShapeCIC2dStatic <PointType>::InCell(x1, x2, gridData->r, gridData->z,
index, CICArray[index])) return searchIndexes[i];
                        }
                };
                break;
        }
        return -1;
};*/
template std::vector<unsigned int> ParticleGridInterface2d<float>::CheckParticlesBoundaries<Particles3dcil<float>*>(
    const std::vector<BoundaryContainer2d<float>>& boundaries, std::vector<ElectrodeCurrent<float>>& conductorList,
    Particles3dcil<float>* state1, Particles3dcil<float>* state2, float dt, float charge, float mass);

template std::vector<unsigned int> ParticleGridInterface2d<float>::CheckParticlesBoundaries<Particles2d<float>*>(
    const std::vector<BoundaryContainer2d<float>>& boundaries, std::vector<ElectrodeCurrent<float>>& conductorList,
    Particles2d<float>* state1, Particles2d<float>* state2, float dt, float charge, float mass);

template std::vector<unsigned int> ParticleGridInterface2d<float>::CheckParticlesBoundaries<Particles2dpolar<float>*>(
    const std::vector<BoundaryContainer2d<float>>& boundaries, std::vector<ElectrodeCurrent<float>>& conductorList,
    Particles2dpolar<float>* state1, Particles2dpolar<float>* state2, float dt, float charge, float mass);

template std::vector<unsigned int> ParticleGridInterface2d<double>::CheckParticlesBoundaries<Particles3dcil<double>*>(
    const std::vector<BoundaryContainer2d<double>>& boundaries, std::vector<ElectrodeCurrent<double>>& conductorList,
    Particles3dcil<double>* state1, Particles3dcil<double>* state2, double dt, double charge, double mass);

template std::vector<unsigned int> ParticleGridInterface2d<double>::CheckParticlesBoundaries<Particles2d<double>*>(
    const std::vector<BoundaryContainer2d<double>>& boundaries, std::vector<ElectrodeCurrent<double>>& conductorList,
    Particles2d<double>* state1, Particles2d<double>* state2, double dt, double charge, double mass);

template std::vector<unsigned int> ParticleGridInterface2d<double>::CheckParticlesBoundaries<Particles2dpolar<double>*>(
    const std::vector<BoundaryContainer2d<double>>& boundaries, std::vector<ElectrodeCurrent<double>>& conductorList,
    Particles2dpolar<double>* state1, Particles2dpolar<double>* state2, double dt, double charge, double mass);

template <class PointType>
template <class particlesType>
std::vector<unsigned int> ParticleGridInterface2d<PointType>::CheckParticlesBoundaries(
    const std::vector<BoundaryContainer2d<PointType>>& boundaries,
    std::vector<ElectrodeCurrent<PointType>>& conductorList, particlesType state1, particlesType state2, PointType dt,
    PointType charge, PointType mass)
{
    int                    tmp;
    DGeo::Point<PointType> tmpPoint;
    int                    flagContinue;

    std::vector<unsigned int>                        result;
    std::vector<std::vector<unsigned int>>           particlesNumbers(s + 1);
    std::vector<std::vector<DGeo::Edge<PointType>>>  intersectionEdges(s + 1);
    std::vector<std::vector<DGeo::Point<PointType>>> intersectionPoints(s + 1);
    std::vector<std::vector<int>>                    intersectionBoundaryNumber(s + 1);

    std::vector<DGeo::Point<PointType>> intersectionPointsConductor;
    std::vector<unsigned int>           particlesNumbersConductor;
    std::vector<int>                    intersectionBoundaryNumberConductor;

    int                   flagNotBoundary;
    int                   index;
    bool                  flagIntersection;
    DGeo::Edge<PointType> traceEdge;

    std::vector<unsigned int> isEmitted;

    for (int k = 0; k < state2->NParticles; k++)
    {

        flagIntersection = false;
        index            = templNumb.getElem(state2->cellsNumbers[k], 0, 0);

        switch (priorityParticleShapeType)
        {
        case 1:
            flagNotBoundary = TSCcellArray.isBoundary[index];
            break;
        case 0:
            flagNotBoundary = CICcellArray.isBoundary[index];
            break;
        }

        if (int(flagNotBoundary) == 0)
            continue;

        flagContinue       = 0;
        traceEdge.point1.x = state1->GetPointerToCartesianX()[k];
        traceEdge.point1.y = state1->GetPointerToCartesianY()[k];
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
                if (int(state2->flagEmitted[k]) != 1)
                {
                    /*flagIntersection = boundaries[list[i][j]].IsIntersection(traceEdge, tmpPoint, tmp);
                    if (flagIntersection)
                    {
                            state2->flagEmitted[k] = 2;
                            flagContinue = 1;
                            break;
                    }*/
                    state2->flagEmitted[k] = 1;
                    flagContinue           = 1;
                    flagIntersection       = false;
                    break;
                }
                else
                {
                    flagIntersection = boundaries[list[i][j]].IsIntersectionLight(traceEdge, tmpPoint, tmp);
                }
                if (flagIntersection)
                {
                    particlesNumbers[i].push_back(k);
                    intersectionEdges[i].push_back(boundaries[list[i][j]].EdgesData[tmp]);
                    intersectionPoints[i].push_back(tmpPoint);
                    intersectionBoundaryNumber[i].push_back(-1);
                    for (int k = 0; k < conductorList.size(); k++)
                    {
                        for (int k1 = 0; k1 < conductorList[k].boundaryNumbers.size(); k1++)
                        {
                            if (list[i][j] == conductorList[k].boundaryNumbers[k1])
                            {
                                intersectionBoundaryNumber[i].back() = k;
                                break;
                            }
                        }
                    }
                    flagContinue = 1;
                    break;
                };
            }
            if (flagContinue)
                break;
        };

        if (flagContinue)
            continue;
    }

    for (int i = 0; i < s; i++)
    {
        ApplyBoundaryCondition(result, boundaryConditions->GetConditionPropertyType(i), conductorList,
                               boundaryConditions->GetConditionProperties(i), state1, state2, particlesNumbers[i],
                               intersectionEdges[i], intersectionPoints[i], intersectionBoundaryNumber[i], dt, charge,
                               mass);
    }

    for (int k = 0; k < state2->NParticles; k++)
    {

        flagContinue = 0;

        flagIntersection = false;
        index            = templNumb.getElem(state2->cellsNumbers[k], 0, 0);

        switch (priorityParticleShapeType)
        {
        case 1:
            flagNotBoundary = TSCcellArray.isBoundary[index];
            break;
        case 0:
            flagNotBoundary = CICcellArray.isBoundary[index];
            break;
        }

        if (int(flagNotBoundary) == 0)
            continue;

        traceEdge.point1.x = state1->GetPointerToCartesianX()[k];
        traceEdge.point1.y = state1->GetPointerToCartesianY()[k];
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

            if (int(state2->flagEmitted[k]) != 1)
            {
                /*flagIntersection = boundaries[list1[j]].IsIntersection(traceEdge, tmpPoint, tmp);
                if (flagIntersection && k == 22)
                {
                        int tt = 0;
                }
                if (flagIntersection)
                {
                        state2->flagEmitted[k] = 1;
                        flagContinue = 1;
                        break;
                }*/

                state2->flagEmitted[k] = 1;
                flagContinue           = 1;
                flagIntersection       = false;
                break;
            }
            else
            {
                flagIntersection = boundaries[list1[j]].IsIntersectionLight(traceEdge, tmpPoint, tmp);
            }

            if (flagIntersection)
            {
                state2->SetCartesianPosition(k, tmpPoint.x, tmpPoint.y);

                particlesNumbers[s].push_back(k);
                intersectionEdges[s].push_back(boundaries[list1[j]].EdgesData[tmp]);
                break;
            };
        }
        if (flagIntersection == false)
            continue;
        flagContinue = 0;
        for (int j = 0; j < conductorList.size(); j++)
        {
            for (int k1 = 0; k1 < conductorList[j].boundaryNumbers.size(); k1++)
            {
                if (boundaries[conductorList[j].boundaryNumbers[k1]].IsIntersectionLight(traceEdge, tmpPoint, tmp))
                {
                    intersectionPointsConductor.push_back(tmpPoint);
                    particlesNumbersConductor.push_back(k);
                    intersectionBoundaryNumberConductor.push_back(j);
                    flagContinue = 1;
                    break;
                }
            };
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
        };
        if (flag == 0)
            result.push_back(particlesNumbers[s][j]);
    }

    //	ApplyDefaultCondition(conductorList, state1, state2, particlesNumbers[s], intersectionEdges[s],
    // intersectionPoints[s], intersectionBoundaryNumber[s]);
    ApplyDefaultCondition(conductorList, state1, state2, particlesNumbersConductor, intersectionEdges[s],
                          intersectionPointsConductor, intersectionBoundaryNumberConductor, dt, charge, mass);

    for (int k = 0; k < state2->NParticles; k++)
    {
        if (int(state2->flagEmitted[k]) == 2)
        {
            state2->flagEmitted[k] = 1;
        }
    }
    return result;
};

template void ParticleGridInterface2d<float>::ApplyBoundaryCondition<Particles3dcil<float>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<float>>& conductorList, std::vector<double> conditionProperties,
    Particles3dcil<float>* state1, Particles3dcil<float>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<float>>& intersectionEdges, const std::vector<DGeo::Point<float>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass);

template void ParticleGridInterface2d<float>::ApplyBoundaryCondition<Particles2d<float>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<float>>& conductorList, std::vector<double> conditionProperties,
    Particles2d<float>* state1, Particles2d<float>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<float>>& intersectionEdges, const std::vector<DGeo::Point<float>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass);

template void ParticleGridInterface2d<float>::ApplyBoundaryCondition<Particles2dpolar<float>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<float>>& conductorList, std::vector<double> conditionProperties,
    Particles2dpolar<float>* state1, Particles2dpolar<float>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<float>>& intersectionEdges, const std::vector<DGeo::Point<float>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass);

template void ParticleGridInterface2d<double>::ApplyBoundaryCondition<Particles3dcil<double>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<double>>& conductorList, std::vector<double> conditionProperties,
    Particles3dcil<double>* state1, Particles3dcil<double>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<double>>&  intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template void ParticleGridInterface2d<double>::ApplyBoundaryCondition<Particles2d<double>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<double>>& conductorList, std::vector<double> conditionProperties,
    Particles2d<double>* state1, Particles2d<double>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<double>>&  intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template void ParticleGridInterface2d<double>::ApplyBoundaryCondition<Particles2dpolar<double>*>(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<double>>& conductorList, std::vector<double> conditionProperties,
    Particles2dpolar<double>* state1, Particles2dpolar<double>* state2,
    const std::vector<unsigned int>& particlesNumbers, const std::vector<DGeo::Edge<double>>& intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::ApplyBoundaryCondition(
    std::vector<unsigned int>& empty, const std::string& conditionType,
    std::vector<ElectrodeCurrent<PointType>>& conductorList, std::vector<double> conditionProperties,
    particlesType state1, particlesType state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
    const std::vector<DGeo::Point<PointType>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    PointType dt, PointType charge, PointType mass)
{
    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[0]) //"Metal Absorbtion"
    {
        for (int i = 0; i < particlesNumbers.size(); i++)
        {
            if (intersectionBoundaryNumber[i] != -1)
                conductorList[intersectionBoundaryNumber[i]].AddCharge(
                    state2->q[particlesNumbers[i]] * (1 - conditionProperties[0]), intersectionPoints[i], dt,
                    state2->GetEnergy(particlesNumbers[i], mass), charge);
            state2->q[particlesNumbers[i]] = state2->q[particlesNumbers[i]] * conditionProperties[0];
            state2->currentCoef[particlesNumbers[i]] =
                state2->currentCoef[particlesNumbers[i]] * conditionProperties[0];

            if (std::abs(state2->q[particlesNumbers[i]]) < 0.01 * state2->avCharge)
            {
                empty.push_back(particlesNumbers[i]);
                state2->SetCartesianPosition(particlesNumbers[i], intersectionPoints[i].x, intersectionPoints[i].y);
            }
        };
    };

    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[2]) //"backscattering"
    {
        for (int k = 0; k < particlesNumbers.size(); k++)
        {
            DGeo::Edge<PointType> traceEdge;
            traceEdge.point1.x = state1->GetPointerToCartesianX()[particlesNumbers[k]];
            traceEdge.point1.y = state1->GetPointerToCartesianY()[particlesNumbers[k]];
            traceEdge.point1.z = 0;

            traceEdge.point2.x   = state2->GetPointerToCartesianX()[particlesNumbers[k]];
            traceEdge.point2.y   = state2->GetPointerToCartesianY()[particlesNumbers[k]];
            traceEdge.point2.z   = 0;
            PointType alphaTrace = traceEdge.alpha();
            PointType alphaEdge  = intersectionEdges[k].alpha();
            PointType dalpha     = alphaEdge - alphaTrace;

            DGeo::Point<PointType> newPoint2 = traceEdge.point2.rotate(2 * dalpha, intersectionPoints[k]);
            state2->SetCartesianPosition(particlesNumbers[k], newPoint2.x, newPoint2.y);

            DGeo::Point<PointType> velocityPoint1 = traceEdge.Middle();
            DGeo::Point<PointType> velocityPoint2;
            velocityPoint2.x = velocityPoint1.x + state2->GetCartesianPX(particlesNumbers[k]);
            velocityPoint2.y = velocityPoint1.y + state2->GetCartesianPY(particlesNumbers[k]);
            velocityPoint2.z = 0;

            DGeo::Point<PointType> newvelocityPoint1 = velocityPoint1.rotate(2 * dalpha, intersectionPoints[k]);
            DGeo::Point<PointType> newvelocityPoint2 = velocityPoint2.rotate(2 * dalpha, intersectionPoints[k]);

            state2->SetCartesianMomentum(particlesNumbers[k], newvelocityPoint2.x - newvelocityPoint1.x,
                                         newvelocityPoint2.y - newvelocityPoint1.y);

            if (intersectionBoundaryNumber[k] != -1)
                conductorList[intersectionBoundaryNumber[k]].AddCharge(
                    state2->q[particlesNumbers[k]] * (1 - conditionProperties[0]), intersectionPoints[k], dt,
                    (1 - conditionProperties[1]) * state2->GetEnergy(particlesNumbers[k], mass), charge);

            state2->MultiplyMomentum(particlesNumbers[k], conditionProperties[1]);

            state2->q[particlesNumbers[k]] = state2->q[particlesNumbers[k]] * conditionProperties[0];

            state2->currentCoef[particlesNumbers[k]] =
                state2->currentCoef[particlesNumbers[k]] * conditionProperties[0];

            if (std::abs(state2->q[particlesNumbers[k]]) < 0.01 * state2->avCharge)
            {
                empty.push_back(particlesNumbers[k]);
                state2->SetCartesianPosition(particlesNumbers[k], intersectionPoints[k].x, intersectionPoints[k].y);
            }
        }
    };

    if (conditionType == flagStringsSolver::flowBoundaryTypeNames[1]) // "Reflection"
    {
        for (int k = 0; k < particlesNumbers.size(); k++)
        {
            DGeo::Edge<PointType> traceEdge;
            traceEdge.point1.x = state1->GetPointerToCartesianX()[particlesNumbers[k]];
            traceEdge.point1.y = state1->GetPointerToCartesianY()[particlesNumbers[k]];
            traceEdge.point1.z = 0;

            traceEdge.point2.x   = state2->GetPointerToCartesianX()[particlesNumbers[k]];
            traceEdge.point2.y   = state2->GetPointerToCartesianY()[particlesNumbers[k]];
            traceEdge.point2.z   = 0;
            PointType alphaTrace = traceEdge.alpha();
            PointType alphaEdge  = intersectionEdges[k].alpha();
            PointType dalpha     = alphaEdge - alphaTrace;

            DGeo::Point<PointType> newPoint2 = traceEdge.point2.rotate(2 * dalpha, intersectionPoints[k]);
            state2->SetCartesianPosition(particlesNumbers[k], newPoint2.x, newPoint2.y);

            DGeo::Point<PointType> velocityPoint1 = traceEdge.Middle();
            DGeo::Point<PointType> velocityPoint2;
            velocityPoint2.x = velocityPoint1.x + state2->GetCartesianPX(particlesNumbers[k]);
            velocityPoint2.y = velocityPoint1.y + state2->GetCartesianPY(particlesNumbers[k]);
            velocityPoint2.z = 0;

            DGeo::Point<PointType> newvelocityPoint1 = velocityPoint1.rotate(2 * dalpha, intersectionPoints[k]);
            DGeo::Point<PointType> newvelocityPoint2 = velocityPoint2.rotate(2 * dalpha, intersectionPoints[k]);

            state2->SetCartesianMomentum(particlesNumbers[k], newvelocityPoint2.x - newvelocityPoint1.x,
                                         newvelocityPoint2.y - newvelocityPoint1.y);
        }
    };
};

template void ParticleGridInterface2d<float>::ApplyDefaultCondition<Particles3dcil<float>*>(
    std::vector<ElectrodeCurrent<float>>& conductorList, Particles3dcil<float>* state1, Particles3dcil<float>* state2,
    const std::vector<unsigned int>& particlesNumbers, const std::vector<DGeo::Edge<float>>& intersectionEdges,
    const std::vector<DGeo::Point<float>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    float dt, float charge, float mass);

template void ParticleGridInterface2d<float>::ApplyDefaultCondition<Particles2d<float>*>(
    std::vector<ElectrodeCurrent<float>>& conductorList, Particles2d<float>* state1, Particles2d<float>* state2,
    const std::vector<unsigned int>& particlesNumbers, const std::vector<DGeo::Edge<float>>& intersectionEdges,
    const std::vector<DGeo::Point<float>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    float dt, float charge, float mass);

template void ParticleGridInterface2d<float>::ApplyDefaultCondition<Particles2dpolar<float>*>(
    std::vector<ElectrodeCurrent<float>>& conductorList, Particles2dpolar<float>* state1,
    Particles2dpolar<float>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<float>>& intersectionEdges, const std::vector<DGeo::Point<float>>& intersectionPoints,
    const std::vector<int>& intersectionBoundaryNumber, float dt, float charge, float mass);

template void ParticleGridInterface2d<double>::ApplyDefaultCondition<Particles3dcil<double>*>(
    std::vector<ElectrodeCurrent<double>>& conductorList, Particles3dcil<double>* state1,
    Particles3dcil<double>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<double>>&  intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template void ParticleGridInterface2d<double>::ApplyDefaultCondition<Particles2d<double>*>(
    std::vector<ElectrodeCurrent<double>>& conductorList, Particles2d<double>* state1, Particles2d<double>* state2,
    const std::vector<unsigned int>& particlesNumbers, const std::vector<DGeo::Edge<double>>& intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template void ParticleGridInterface2d<double>::ApplyDefaultCondition<Particles2dpolar<double>*>(
    std::vector<ElectrodeCurrent<double>>& conductorList, Particles2dpolar<double>* state1,
    Particles2dpolar<double>* state2, const std::vector<unsigned int>& particlesNumbers,
    const std::vector<DGeo::Edge<double>>&  intersectionEdges,
    const std::vector<DGeo::Point<double>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    double dt, double charge, double mass);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::ApplyDefaultCondition(
    std::vector<ElectrodeCurrent<PointType>>& conductorList, particlesType state1, particlesType state2,
    const std::vector<unsigned int>& particlesNumbers, const std::vector<DGeo::Edge<PointType>>& intersectionEdges,
    const std::vector<DGeo::Point<PointType>>& intersectionPoints, const std::vector<int>& intersectionBoundaryNumber,
    PointType dt, PointType charge, PointType mass)
{
    for (int i = 0; i < particlesNumbers.size(); i++)
    {
        if (intersectionBoundaryNumber[i] != -1)
            conductorList[intersectionBoundaryNumber[i]].AddCharge(
                state2->q[particlesNumbers[i]], intersectionPoints[i], dt, state2->GetEnergy(particlesNumbers[i], mass),
                charge);
    }
    //	state2->EmptyPlaces = particlesNumbers;
};

template <class PointType>
void ParticleGridInterface2d<PointType>::Grid2Particles(Particles2d<PointType>*            particles,
                                                        const GridData2d<PointType>* const gridData, PointType (*W)[9])
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray.ValueInterpolate(index, &W[i][0], gridData->Ex, particles->Ex[i]);
            TSCcellArray.ValueInterpolate(index, &W[i][0], gridData->Ey, particles->Ey[i]);
        }
        break;
    case 0:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            CICcellArray.ValueInterpolate(index, &W[i][0], gridData->Ex, particles->Ex[i]);
            CICcellArray.ValueInterpolate(index, &W[i][0], gridData->Ey, particles->Ey[i]);
            particles->cellSize[i] = CICcellArray.GetSize(index);
        }
    }
};
template <class PointType>
void ParticleGridInterface2d<PointType>::Grid2Particles(Particles3dcil<PointType>*            particles,
                                                        const GridData2daxs<PointType>* const gridData,
                                                        PointType (*W)[9])
{
    int index;
    switch (priorityParticleShapeType)
    {
    case 1:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray.ValueInterpolate(index, &W[i][0], gridData->Er, particles->Er[i]);
            TSCcellArray.ValueInterpolate(index, &W[i][0], gridData->Ez, particles->Ez[i]);
            TSCcellArray.ValueInterpolate(index, &W[i][0], gridData->Bphi, particles->Bphi[i]);
        }
        break;
    case 0:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            CICcellArray.ValueInterpolate(index, &W[i][0], gridData->Er, particles->Er[i]);
            CICcellArray.ValueInterpolate(index, &W[i][0], gridData->Ez, particles->Ez[i]);
            CICcellArray.ValueInterpolate(index, &W[i][0], gridData->Bphi, particles->Bphi[i]);
            particles->cellSize[i] = CICcellArray.GetSize(index);
        }
    }
};
template <class PointType>
void ParticleGridInterface2d<PointType>::Particles2Grid(const Particles2d<PointType>* const particles,
                                                        GridData2d<PointType>* gridData, PointType (*W)[9]){

};
template <class PointType>
void ParticleGridInterface2d<PointType>::Particles2Grid(const Particles3dcil<PointType>* const particles,
                                                        GridData2daxs<PointType>* gridData, PointType (*W)[9])
{
    int index;

    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            TSCcellArray.ChargeCalculate(index, &W[i][0], particles->q[i], gridData->rho);
        }
        break;

    case 0:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            CICcellArray.ChargeCalculate(index, &W[i][0], particles->q[i], gridData->rho);
        }
        break;
    }
};

template void ParticleGridInterface2d<float>::Particles2GridPTI<Particles3dcil<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles3dcil<float>* state1, float (*W1)[9],
    Particles3dcil<float>* state2, float (*W2)[9], float* rho, float dt);

template void ParticleGridInterface2d<float>::Particles2GridPTI<Particles2d<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles2d<float>* state1, float (*W1)[9],
    Particles2d<float>* state2, float (*W2)[9], float* rho, float dt);

template void ParticleGridInterface2d<float>::Particles2GridPTI<Particles2dpolar<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles2dpolar<float>* state1, float (*W1)[9],
    Particles2dpolar<float>* state2, float (*W2)[9], float* rho, float dt);

template void ParticleGridInterface2d<double>::Particles2GridPTI<Particles3dcil<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles3dcil<double>* state1, double (*W1)[9],
    Particles3dcil<double>* state2, double (*W2)[9], double* rho, double dt);

template void ParticleGridInterface2d<double>::Particles2GridPTI<Particles2d<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles2d<double>* state1, double (*W1)[9],
    Particles2d<double>* state2, double (*W2)[9], double* rho, double dt);

template void ParticleGridInterface2d<double>::Particles2GridPTI<Particles2dpolar<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, Particles2dpolar<double>* state1, double (*W1)[9],
    Particles2dpolar<double>* state2, double (*W2)[9], double* rho, double dt);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::Particles2GridPTI(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType, arma::mat& Icoef,
    const std::vector<unsigned int>& emissionCells, particlesType state1, PointType (*W1)[9], particlesType state2,
    PointType (*W2)[9], PointType* rho, PointType dt)
{
    int index1;
    int index2;

    PointType* Xar1 = state1->GetPointerToPosition1();
    PointType* Yar1 = state1->GetPointerToPosition2();

    PointType* Xar2 = state2->GetPointerToPosition1();
    PointType* Yar2 = state2->GetPointerToPosition2();

    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = 0; i < state1->NParticles; i++)
        {
            index1 = templNumb.getElem(state1->cellsNumbers[i], 0, 0);
            index2 = templNumb.getElem(state2->cellsNumbers[i], 0, 0);

            TSCcellArray.ChargeCalculate(Xar1[i], Yar1[i], index1, &W1[i][0], Xar2[i], Yar2[i], index2, &W1[i][0], 0, 0,
                                         state2->q[i], dt, rho);
        }
        break;

    case 0:
        for (int i = 0; i < state1->NParticles; i++)
        {
            index1 = templNumb.getElem(state1->cellsNumbers[i], 0, 0);
            index2 = templNumb.getElem(state2->cellsNumbers[i], 0, 0);

            CICcellArray.ChargeCalculate(nearCathodeVolumes, emType, Icoef, state2->currentCoef[i],
                                         state2->flagEmitted[i], emissionCells, state2->startCellNumbersGrid[i],
                                         Xar1[i], Yar1[i], index1, &W1[i][0], Xar2[i], Yar2[i], index2, &W1[i][0], 0, 0,
                                         state2->q[i], dt, rho);
        }
        break;
    }
};

template void ParticleGridInterface2d<float>::InitEmCells<Particles3dcil<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
    Particles3dcil<float>* state1, double dH);
template void ParticleGridInterface2d<float>::InitEmCells<Particles2d<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
    Particles2d<float>* state1, double dH);
template void ParticleGridInterface2d<float>::InitEmCells<Particles2dpolar<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, int emType,
    Particles2dpolar<float>* state1, double dH);

template void ParticleGridInterface2d<double>::InitEmCells<Particles3dcil<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
    Particles3dcil<double>* state1, double dH);
template void ParticleGridInterface2d<double>::InitEmCells<Particles2d<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
    Particles2d<double>* state1, double dH);
template void ParticleGridInterface2d<double>::InitEmCells<Particles2dpolar<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, int emType,
    Particles2dpolar<double>* state1, double dH);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::InitEmCells(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType, particlesType state1,
    double dH)
{
    PointType* Xar = state1->GetPointerToCartesianX();
    PointType* Yar = state1->GetPointerToCartesianY();

    int nrow = templNumb.GetNrow();

    std::vector<int> searchIndexesL(9);

    for (int i = 0; i < state1->NParticles; i++)
    {
        state1->startCellNumbersGrid[i] = -1;
        double px                       = state1->GetCartesianPX(i);
        double py                       = state1->GetCartesianPY(i);
        double p                        = sqrt(px * px + py * py);
        double x1                       = Xar[i];
        double y1                       = Yar[i];

        x1 = x1 + dH * 0.02 * px / p;
        y1 = y1 + dH * 0.02 * py / p;

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

            int index = InCellWithEps(x1, y1, searchIndexesL);

            if (index != -1)
            {
                int index1 = templNumb.getElem(index, 0, 0);

                int n = CICcellArray.InitEmCells(nearCathodeVolumes, x1, y1, index1);
                if (state1->startCellNumbersGrid[i] == -1)
                    state1->startCellNumbersGrid[i] = n;
            }
            x1 = x1 + dH * 0.1 * px / p;
            y1 = y1 + dH * 0.1 * py / p;
        }
        if (state1->startCellNumbersGrid[i] == -1)
            state1->startCellNumbersGrid[i] = state1->startCellNumbersGrid[i - 1];
    }
};

template <class PointType> void ParticleGridInterface2d<PointType>::Charge2Density(PointType* rho)
{
    int index;

    switch (priorityParticleShapeType)
    {
    case 1:

        for (int j = 0; j < TSCcellArray.cellVolume.size(); j++)
            rho[j] = rho[j] / TSCcellArray.CellVolume(j);
        break;

    case 0:
        for (int j = 0; j < CICcellArray.cellVolume.size(); j++)
            rho[j] = rho[j] / CICcellArray.CellVolume(j);
        break;
    }
};

template <class PointType>
void ParticleGridInterface2d<PointType>::Wcalculate(Particles3dcil<PointType>* particles, PointType (*W)[9])
{
    int index;
    indexes.resize(0);
    cells.resize(0);
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            for (int k  = 0; k < 9; k++)
                W[i][k] = 0;
            /*if (TSCcellArray.cellType[index] == 0)
            {
                    indexes.push_back(i);
                    cells.push_back(index);
            }
            else*/
            TSCcellArray.Wcalculate(index, particles->r[i], particles->z[i], &W[i][0]);
        };
        //	WcalculateVector1(templNumb, &TSCcellArray.cellType[0], &particles->cellsNumbers[0], &particles->r[0],
        //&particles->z[0], W, particles->NParticles, x1Array, x2Array, &TSCcellArray.levelHigh[0]);
        //#pragma noinline
        // if (indexes.size()>0)
        // WcalculateVector(&cells[0],&indexes[0], &particles->r[0], &particles->z[0], W, indexes.size(), x1Array,
        // x2Array, &TSCcellArray.levelHigh[0]);
        break;
    case 0:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            for (int k  = 0; k < 9; k++)
                W[i][k] = 0;
            CICcellArray.Wcalculate(index, particles->r[i], particles->z[i], &W[i][0]);
        }

        break;
    }
};
template <class PointType>
void ParticleGridInterface2d<PointType>::Wcalculate(Particles2d<PointType>* particles, PointType (*W)[9])
{
    int index;
    indexes.resize(0);
    cells.resize(0);
    switch (priorityParticleShapeType)
    {
    case 1:

        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            for (int k  = 0; k < 9; k++)
                W[i][k] = 0;
            TSCcellArray.Wcalculate(index, particles->x[i], particles->y[i], &W[i][0]);
        };
        break;
    case 0:
        for (int i = 0; i < particles->NParticles; i++)
        {
            index = templNumb.getElem(particles->cellsNumbers[i], 0, 0);
            for (int k  = 0; k < 9; k++)
                W[i][k] = 0;
            CICcellArray.Wcalculate(index, particles->x[i], particles->y[i], &W[i][0]);
        }

        break;
    }
};

template void
ParticleGridInterface2d<float>::SearchStartCells<Particles3dcil<float>*>(Particles3dcil<float>* particles);

template void
ParticleGridInterface2d<double>::SearchStartCells<Particles3dcil<double>*>(Particles3dcil<double>* particles);

template void ParticleGridInterface2d<float>::SearchStartCells<Particles2d<float>*>(Particles2d<float>* particles);

template void ParticleGridInterface2d<double>::SearchStartCells<Particles2d<double>*>(Particles2d<double>* particles);

template void
ParticleGridInterface2d<float>::SearchStartCells<Particles2dpolar<float>*>(Particles2dpolar<float>* particles);

template void
ParticleGridInterface2d<double>::SearchStartCells<Particles2dpolar<double>*>(Particles2dpolar<double>* particles);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::SearchStartCells(particlesType particles)
{
    particles->startCellNumbers.clear();
    std::vector<int> searchIndexesLoc;

    searchIndexesLoc.resize(templNumb.ArSize());
    for (int k              = 0; k < templNumb.ArSize(); k++)
        searchIndexesLoc[k] = k;

    particles->cellsNumbers.clear();
    particles->cellsNumbers.resize(particles->NParticles);

    particles->startCellNumbers.clear();
    particles->startCellNumbers.resize(particles->NParticles);

    for (int i = 0; i < particles->NParticles; i++)
    {
        particles->cellsNumbers[i] = InCellWithEps(particles->GetPointerToPosition1()[i],
                                                   particles->GetPointerToPosition2()[i], searchIndexesLoc);

        if (particles->flagEmitted[i] == 0)
        {
            particles->startCellNumbers[i] = InCellWithEps(particles->GetPointerToPosition1()[i],
                                                           particles->GetPointerToPosition2()[i], searchIndexesLoc);
            particles->cellsNumbers[i] = particles->startCellNumbers[i];
        };
    }

    /*std::vector<int>  searchIndexesParticles (particles->startCellNumbers.size()*9);

    int nrow = templNumb.GetNrow();
    int base;

    for (int j = 0; j < particles->startCellNumbers.size(); j++)
    {
            base = particles->startCellNumbers[j];
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
    particles->searchIndexesAll = searchIndexesParticles;*/
}

template void ParticleGridInterface2d<float>::SearchStartCellsEmission<Particles3dcil<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, Particles3dcil<float>* particles);

template void ParticleGridInterface2d<double>::SearchStartCellsEmission<Particles3dcil<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, Particles3dcil<double>* particles);

template void ParticleGridInterface2d<float>::SearchStartCellsEmission<Particles2d<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, Particles2d<float>* particles);

template void ParticleGridInterface2d<double>::SearchStartCellsEmission<Particles2d<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, Particles2d<double>* particles);

template void ParticleGridInterface2d<float>::SearchStartCellsEmission<Particles2dpolar<float>*>(
    const std::vector<std::vector<NearCathodeVolume<float>>>& nearCathodeVolumes, Particles2dpolar<float>* particles);

template void ParticleGridInterface2d<double>::SearchStartCellsEmission<Particles2dpolar<double>*>(
    const std::vector<std::vector<NearCathodeVolume<double>>>& nearCathodeVolumes, Particles2dpolar<double>* particles);

template <class PointType>
template <class particlesType>
void ParticleGridInterface2d<PointType>::SearchStartCellsEmission(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, particlesType particles)
{

    PointType currentFrom_dl = 0;
    int       i0             = 0;
    int       i1             = -1;

    while (1)
    {
        currentFrom_dl = 0;
        do
        {
            i1++;
            currentFrom_dl = currentFrom_dl + particles->q[i1];
            if (i1 == particles->NParticles - 1)
                break;
        } while (particles->startCellNumbersGrid[i1] == particles->startCellNumbersGrid[i1 + 1]);

        for (int k = i0; k < i1 + 1; k++)
        {
            particles->currentCoef[k] = particles->q[k] / std::abs(currentFrom_dl);
            if (std::abs(currentFrom_dl) < 1e-17)
                particles->currentCoef[k] = 0;
        }

        if (i1 == particles->NParticles - 1)
            break;

        i0 = i1 + 1;
        s++;
    };
};

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCell<Particles3dcil<float>*>(Particles3dcil<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCell<Particles3dcil<double>*>(Particles3dcil<double>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCell<Particles2d<float>*>(Particles2d<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCell<Particles2d<double>*>(Particles2d<double>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCell<Particles2dpolar<float>*>(Particles2dpolar<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCell<Particles2dpolar<double>*>(Particles2dpolar<double>* particles);

template <class PointType>
template <class particlesType>
std::vector<unsigned int> ParticleGridInterface2d<PointType>::InCell(particlesType particles)
{
    std::vector<unsigned int> remove;
    int                       nrow = templNumb.GetNrow();

    int base;
    int k = 0;

    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();
    ;

    for (int i = 0; i < particles->NParticles; i++)
    {
        if (particles->cellsNumbers[i] == 0)
        {
            particles->cellsNumbers[i] = InCell(Xar[i], Yar[i], particles->startCellNumbers);
        }
        else
        {
            int base         = particles->cellsNumbers[i];
            searchIndexes[0] = base;
            searchIndexes[1] = base - 1;
            searchIndexes[2] = base + 1;

            searchIndexes[3] = base + nrow;
            searchIndexes[4] = base - 1 + nrow;
            searchIndexes[5] = base + 1 + nrow;

            searchIndexes[6] = base - nrow;
            searchIndexes[7] = base - 1 - nrow;
            searchIndexes[8] = base + 1 - nrow;

            particles->cellsNumbers[i] = InCell(Xar[i], Yar[i], searchIndexes);
        };

        if (particles->cellsNumbers[i] == -1)
            remove.push_back(i);

        if (particles->cellsNumbers[i] == -1)
        {
            int a = 0;
        };
    };
    //	particles->removeParticle(remove);
    iter++;
    return remove;
};

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCellWithEps<Particles3dcil<float>*>(Particles3dcil<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCellWithEps<Particles3dcil<double>*>(Particles3dcil<double>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCellWithEps<Particles2d<float>*>(Particles2d<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCellWithEps<Particles2d<double>*>(Particles2d<double>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<float>::InCellWithEps<Particles2dpolar<float>*>(Particles2dpolar<float>* particles);

template std::vector<unsigned int>
ParticleGridInterface2d<double>::InCellWithEps<Particles2dpolar<double>*>(Particles2dpolar<double>* particles);

template <class PointType>
template <class particlesType>
std::vector<unsigned int> ParticleGridInterface2d<PointType>::InCellWithEps(particlesType particles)
{
    std::vector<unsigned int> remove;
    int                       nrow = templNumb.GetNrow();

    int base;
    int k = 0;

    PointType* Xar = particles->GetPointerToPosition1();
    PointType* Yar = particles->GetPointerToPosition2();
    ;

    for (int i = 0; i < particles->NParticles; i++)
    {
        if (particles->cellsNumbers[i] == 0)
        {
            k++;
            particles->cellsNumbers[i] = InCellWithEps(Xar[i], Yar[i], particles->startCellNumbers);
        }
        else
        {
            if (int(particles->flagEmitted[k]) != 1)
                particles->flagEmitted[k] = 1;

            int base         = particles->cellsNumbers[i];
            searchIndexes[0] = base;
            searchIndexes[1] = base - 1;
            searchIndexes[2] = base + 1;

            searchIndexes[3] = base + nrow;
            searchIndexes[4] = base - 1 + nrow;
            searchIndexes[5] = base + 1 + nrow;

            searchIndexes[6] = base - nrow;
            searchIndexes[7] = base - 1 - nrow;
            searchIndexes[8] = base + 1 - nrow;

            particles->cellsNumbers[i] = InCellWithEps(Xar[i], Yar[i], searchIndexes);
        };

        if (particles->cellsNumbers[i] == -1)
            remove.push_back(i);

        if (particles->cellsNumbers[i] == -1)
        {
            int a = 0;
        };
    };
    //	particles->removeParticle(remove);
    iter++;
    return remove;
};
