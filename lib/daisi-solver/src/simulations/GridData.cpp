#include "GridData.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Geom.h"


template <class PointType>
template <class Archive>
void GridDataBase<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& VA;
    ar& EA;
    ar& rho[0];
    ar& X;
    ar& B;
    ar& flagOut;
    ar& CICArray;
}

template <class PointType>
template <class Archive>
void GridDataBase<PointType>::load(Archive& ar, const unsigned int)
{
    rho.resize(1);
    ar& VA;
    ar& EA;
    ar& rho[0];
    ar& X;
    ar& B;
    ar& flagOut;
    ar& CICArray;
    Init();
}

template <class PointType>
void GridDataBase<PointType>::Init()
{
    V = VA;
    VCharge.resize(VA.size());
    E = EA;
    ECol.resize(EA.size());
    for (int i = 0; i < ECol.size(); i++)
        ECol[i].resize(EA[i].size());
}

template <class PointType>
void GridDataBase<PointType>::GetDataIntFlag(void* Array[1], int& size, int& sizeElement, int flag,
                                             int PlotTypeFlag) const
{
    size        = int(E[0].size());
    sizeElement = sizeof(E[0][0]);

    if (PlotTypeFlag == 0)
    {
        if (flag == 0)
            Array[0] = (void*)(&E[0][0]);

        if (flag == 1)
            Array[0] = (void*)(&E[1][0]);
    }

    if (PlotTypeFlag == 1)
    {
        if (flag == 0)
            Array[0] = (void*)(&ECol[0][0]);

        if (flag == 1)
            Array[0] = (void*)(&ECol[1][0]);
    }

	if(PlotTypeFlag == 3)
	{
		if (flag == 0)
			Array[0] = (void*)(&X[0][0]);

		if (flag == 1)
			Array[0] = (void*)(&X[1][0]);
	}
}

template <class PointType>
void GridDataBase<PointType>::ZeroingFieldsBase()
{
    if (VCharge.size() != 0)
    {
        memset(&VCharge[0], 0, VCharge.size() * sizeof(VCharge[0]));
        memset(&VA[0], 0, VA.size() * sizeof(VA[0]));
        memset(&V[0], 0, V.size() * sizeof(V[0]));
    }
}

template <class PointType>
GridDataBase<PointType>::GridDataBase()
{
    rho.resize(1);
}

template <class PointType>
void GridDataBase<PointType>::ZeroingFields(){

}

template <class PointType>
void GridDataBase<PointType>::ApplyTimeDepending(PointType frequency, PointType phase,
                                                 PointType time, std::vector<int>& nonZeros)
{

    PointType k = std::cos(2 * PI() * frequency * time + phase);
}

template <class PointType>
void GridDataBase<PointType>::ApplyTimeDepending(const std::vector<double>& globalPar, double time)
{
    PointType k;
    if (globalPar[0] < 1e-7)
        k = std::cos(globalPar[1]);
    else
        k = std::cos(2 * PI() * globalPar[0] * time + globalPar[1]);

    for (int j = 0; j < E.size(); j++)
        for (int i  = 0; i < E[j].size(); i++)
            E[j][i] = EA[j][i] * k;

    //+ ECol[j][i];

    for (int i = 0; i < V.size(); i++)
        V[i]   = VA[i] * k;

    //+ VCharge[i];
}

template <class PointType>
void GridDataBase<PointType>::ApplyTimeDepending(PointType frequency, PointType phase,
                                                 PointType time){
    /*if (frequency < 1e-3)
    return;

    PointType k = std::cos(2 * PI()*frequency*time + phase);
    for (int i = 0; i < E[0].size(); i++)
    {
    E[0][i] = ErA[i] * k + ErCol[i];
    Ez[i] = EzA[i] * k + EzCol[i];
    //V[i] = VA[i] + VCharge[i];
    };*/

}

template <class PointType>
void GridDataBase<PointType>::Clear()
{
    X.clear();
    E.clear();
    EA.clear();
    ECol.clear();

    B.clear();

    rho.clear();
    V.clear();
    VA.clear();
    VCharge.clear();
    flagOut.clear();
    CICArray.clear();
}

template <class PointType>
void GridDataBase<PointType>::SetCells(Dmath::imat& flagMatrix, Dmath::imat& templNumb,
                                       int nVertexX, int nVertexY, int problemType)
{

    int levelHigh = 0;
    for (int j = 1; j < nVertexY - 1; j++)
    {

        for (int i = 1; i < nVertexX - 1; i++)
        {
            if (templNumb(i, j) != -1)
            {
                flagOut.push_back(flagMatrix(i, j));
                if (templNumb(i + 1, j) != -1 && templNumb(i + 1, j + 1) != -1 &&
                    templNumb(i, j + 1) != -1)
                {
                    levelHigh = templNumb(i, j + 1);
                    CICArray.push_back(levelHigh);
                }
                else
                    CICArray.push_back(-1);
            }
        }
    }
}

template <class PointType>
void GridDataBase<PointType>::Init(std::vector<DGeo::Point<PointType>> serialMeshData,
                                   int                                 problemType)
{
    Clear();
    rho.resize(1);
    switch (problemType)
    {
    case 1:
        X.resize(2);
        E.resize(2);
        EA.resize(2);
        ECol.resize(2);
        break;
    case 2:
        X.resize(2);
        E.resize(2);
        EA.resize(2);
        ECol.resize(2);
        B.resize(1);
        break;
    case 3:
        X.resize(2);
        E.resize(2);
        EA.resize(2);
        ECol.resize(2);
        break;
    }
    int n = serialMeshData.size();
    for (int i = 0; i < X.size(); i++)
    {
        X[i].resize(n);
        E[i].resize(n);
        EA[i].resize(n);
        ECol[i].resize(n);
    }
    for (int i = 0; i < B.size(); i++)
        B[i].resize(n);

    V.resize(n);
    VA.resize(n);
    VCharge.resize(n);
    rho[0].resize(n);

    if (problemType == 1 || problemType == 2)
    {

        for (int i = 0; i < serialMeshData.size(); i++)
        {
            X[0][i] = serialMeshData[i].x;
            X[1][i] = serialMeshData[i].y;
        }
    }
    if (problemType == 3)
    {
        for (int i = 0; i < serialMeshData.size(); i++)
            Dmath::Cartesian2Polar(serialMeshData[i].x, serialMeshData[i].y, X[0][i], X[1][i]);
    }
}

template <class PointType>
void GridDataBase<PointType>::Summrho()
{
    for (int i = 1; i < rho.size(); i++)
    {
        for (int j    = 0; j < rho[i].size(); j++)
            rho[0][j] = rho[0][j] + rho[i][j];
        /*if (sizeof(PointType) == 4)
                vsAdd(rho[0].size(), (float*)&rho[0][0], (float*)&rho[i][0], (float*)&rho[0][0]);
        if (sizeof(PointType) == 8)
                vdAdd(rho[0].size(), (double*)&rho[0][0], (double*)&rho[i][0],
        (double*)&rho[0][0]);*/
    }
}

template <class PointType>
void GridDataBase<PointType>::densityReset()
{
    for (int i = 0; i < rho.size(); i++)
    {
        memset(&rho[i][0], 0, rho[i].size() * sizeof(rho[0][0]));
    }
}

template <class PointType>
std::vector<PointType>& GridDataBase<PointType>::Getrho()
{
    return rho[0];
}

template <class PointType>
std::vector<PointType>& GridDataBase<PointType>::GetV()
{
    return V;
}

template <class PointType>
std::vector<PointType>& GridDataBase<PointType>::GetVCharge()
{
    return VCharge;
}

template <class PointType>
std::vector<PointType>& GridDataBase<PointType>::GetVA()
{
    return VA;
}

template <class PointType>
void GridDataBase<PointType>::InitParallel(int numThreads)
{
    rho.resize(numThreads);
    for (int i = 1; i < numThreads; i++)
        rho[i].resize(rho[0].size());
}

template <class PointType>
std::vector<PointType>& GridDataBase<PointType>::Getrho(int thread)
{
    return rho[thread];
}

template <class PointType>
NearCathodeVolume<PointType>::NearCathodeVolume(DGeo::Edge<PointType> Edge1In,
                                                DGeo::Edge<PointType> Edge2In, int flag, int flagEm)
{
    Edge.resize(4);
    // ��� �������� Edge ������ ��������� �������, ��� ��������� �������
    Edge[0]        = Edge1In;
    Edge[1]        = Edge2In;
    Edge[2].point1 = Edge[0].point1;
    Edge[2].point2 = Edge[1].point1;

    Edge[3].point1 = Edge[0].point2;
    Edge[3].point2 = Edge[1].point2;

    flagCurrentLimited = flagEm;

    for (int i = 0; i < 4; i++)
    {
        fieldPointsX[i] = Edge[i].Middle().x;
        fieldPointsY[i] = Edge[i].Middle().y;

        if (flag == 1 || flag == 3)
            Areas[i] = Edge[i].length();
        if (flag == 2)
        {
            if (std::abs(Edge[i].point1.y - Edge[i].point2.y) < 1e-12)
                Areas[i] = PI() * std::abs(Edge[i].point1.x * Edge[i].point1.x -
                                                 Edge[i].point2.x * Edge[i].point2.x);
            else
                Areas[i] =
                    PI() * (Edge[i].point1.x + Edge[i].point2.x) * Edge[i].length();
        }
    }
    /*	if (flag == 2)
            {
                    if (std::abs(Edge[i].point1.y - Edge[i].point2.y) < 1e-12)
                            volume = PI()*std::abs(Edge[i].point1.x*Edge[i].point1.x -
       Edge[i].point2.x* Edge[i].point2.x); else Areas[i] = PI()*(Edge[i].point1.x +
       Edge[i].point2.x)*Edge[i].length();
            }*/

    if (flag == 1 || flag == 3)
        volume = Edge[0].length() * Edge[1].length();

    PointType              h   = Edge[1].length();
    PointType              eps = h * 1e-5;
    DGeo::Point<PointType> ptmp;

    if (Edge[0].IsOnOneSide(Edge[0].GetNormalPoint1(h), Edge[1].Middle(), eps))
        ptmp = Edge[0].GetNormalPoint2(1.0);
    else
        ptmp = Edge[0].GetNormalPoint1(1.0);

    normalX[0] = ptmp.x - Edge[0].Middle().x;
    normalY[0] = ptmp.y - Edge[0].Middle().y;

    if (Edge[1].IsOnOneSide(Edge[1].GetNormalPoint1(h), Edge[0].Middle(), eps))
    {
        ptmp = Edge[1].GetNormalPoint2(1.0);
    }
    else
        ptmp = Edge[1].GetNormalPoint1(1.0);

    normalX[1] = ptmp.x - Edge[1].Middle().x;
    normalY[1] = ptmp.y - Edge[1].Middle().y;

    if (Edge[2].IsOnOneSide(Edge[2].GetNormalPoint1(h), Edge[3].Middle(), eps))
    {
        ptmp = Edge[2].GetNormalPoint2(1.0);
    }
    else
        ptmp = Edge[2].GetNormalPoint1(1.0);

    normalX[2] = ptmp.x - Edge[2].Middle().x;
    normalY[2] = ptmp.y - Edge[2].Middle().y;

    if (Edge[3].IsOnOneSide(Edge[3].GetNormalPoint1(h), Edge[2].Middle(), eps))
    {
        ptmp = Edge[3].GetNormalPoint2(1.0);
    }
    else
        ptmp = Edge[3].GetNormalPoint1(1.0);

    normalX[3] = ptmp.x - Edge[3].Middle().x;
    normalY[3] = ptmp.y - Edge[3].Middle().y;
}

template <class PointType>
bool NearCathodeVolume<PointType>::InCell(PointType x1, PointType x2) const
{

    DGeo::Point<PointType> ptmp;
    ptmp.x = x1;
    ptmp.y = x2;
    ptmp.z = 0;

    PointType h   = Edge[1].length();
    PointType eps = h * 1e-5;

    for (int i = 0; i < 4; i++)
    {
        if (Edge[i].IsOnEdge(ptmp, eps) < 2)
            return true;
    }

    DGeo::Edge<PointType> etmp;
    etmp.point1 = ptmp;

    DGeo::Point<PointType> ptmp1;
    ptmp1.x = x1 + 10 * h;
    ptmp1.y = x2 + 10 * h;
    ptmp1.z = 0;

    etmp.point2 = ptmp1;

    int inters = 0;
    for (int i = 0; i < 4; i++)
    {
        if (Edge[i].IsIntersectionEdgeLight(etmp, eps, &ptmp1))
            inters++;
    }
    if (inters % 2)
        return true;

    return false;
}

template <class PointType>
bool GridDataBase<PointType>::SearchIntersectionWithEdge(int                     cellNumb,
                                                         DGeo::Edge<PointType>   EdgeIntersection,
                                                         DGeo::Point<PointType>* startPoint) const
{
    std::vector<DGeo::Edge<PointType>> edges = GetCellEdgesArray(cellNumb);

    for (int i = 0; i < edges.size(); i++)
    {
        if (edges[i].IsIntersectionEdge(EdgeIntersection, EdgeIntersection.length() * 1e-5,
                                        startPoint) &&
            !edges[i].IsParallel(EdgeIntersection, EdgeIntersection.length() * 1e-5))
            return true;
    }

    return false;
}


template class GridDataBase<double>;
template class GridDataBase<float>;
template class NearCathodeVolume<double>;
template class NearCathodeVolume<float>;

template void
GridDataBase<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);
template void
GridDataBase<double>::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                            const unsigned int file_version) const;

template void
GridDataBase<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                            const unsigned int file_version);

template void
GridDataBase<float>::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version) const;