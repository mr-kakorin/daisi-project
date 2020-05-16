#include "GridData.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "ParticleShape2d.h"

template class GridData3d<double>;
template class GridData3d<float>;

template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Getx()
{
    return this->X[0];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Gety()
{
    return this->X[1];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Getz()
{
    return this->X[2];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_Ex()
{
    return this->E[0];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_Ey()
{
    return this->E[1];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_Ez()
{
    return this->E[2];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_ExA()
{
    return this->EA[0];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_EyA()
{
    return this->EA[1];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_EzA()
{
    return this->EA[2];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_ExCol()
{
    return this->ECol[0];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_EyCol()
{
    return this->ECol[1];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::Get_EzCol()
{
    return this->ECol[2];
};

template <class PointType>
std::vector<PointType>& GridData3d<PointType>::GetBx()
{
    return this->B[0];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::GetBy()
{
    return this->B[1];
};
template <class PointType>
std::vector<PointType>& GridData3d<PointType>::GetBz()
{
    return this->B[2];
};

template <class PointType>
float GridData3d<PointType>::GetMaxSixe() const
{
    float Hmax = 0;
    float h;
    for (int i = 1; i < this->X[0].size(); i++)
    {

        h = std::abs(this->X[1][i] - this->X[1][i - 1]);

        if (h > Hmax)
            Hmax = h;
    }

    return float(std::abs(this->X[0][1] - this->X[0][0]));

    return std::max(float(std::abs(this->X[0][1] - this->X[0][0])), float(Hmax));
};

template <class PointType>
int GridData3d<PointType>::InCell(double x1, double x2, double x3) const
{
    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCell(PointType(x1), PointType(x2), PointType(x3),
                                                        this->X[0], this->X[1], this->X[2], i,
                                                        this->CICArray, this->CICArrayZ))
            break;
    }
    if (i == this->CICArray.size() - 1)
        return -1;
    return i;
};

template <class PointType>
void GridData3d<PointType>::GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                                    int PlotTypeFlag) const {

};

template <class PointType>
float GridData3d<PointType>::interpolatePoint(double x1, double x2, double x3, std::string value,
                                              int PlotTypeFlag) const
{
    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCell(PointType(x1), PointType(x2), PointType(x3),
                                                        this->X[0], this->X[1], this->X[2], i,
                                                        this->CICArray, this->CICArrayZ))
            break;
    }
    if (i == this->CICArray.size() - 1)
        return NAN;

    std::vector<PointType> W(8);
    ParticleShapeCIC2dStatic<PointType>::Wcalculate3d(PointType(x1), PointType(x2), PointType(x3),
                                                      this->X[0], this->X[1], this->X[2], W, i,
                                                      this->CICArray, this->CICArrayZ);

    PointType result = 0;

    if (value == flagStringsSolver::PlotFlags3d[0])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[0], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[1])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[1], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[2])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[2], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[3])
    {
        PointType tmp1;
        PointType tmp2;
        PointType tmp3;

        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[0], tmp1, i,
                                                                this->CICArray, this->CICArrayZ);
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[1], tmp2, i,
                                                                this->CICArray, this->CICArrayZ);
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->E[2], tmp3, i,
                                                                this->CICArray, this->CICArrayZ);

        result = sqrt(tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3);
    }

    if (value == flagStringsSolver::PlotFlags3d[4])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->B[0], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[5])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->B[1], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[6])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->B[2], result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[7])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->V, result, i,
                                                                this->CICArray, this->CICArrayZ);

    if (value == flagStringsSolver::PlotFlags3d[8])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(W, this->rho[0], result, i,
                                                                this->CICArray, this->CICArrayZ);

    return result;
};

template <class PointType>
void GridData3d<PointType>::interpolatePoint(double x1, double x2, double x3, double& Exin,
                                             double& Eyin) const
{
}

template <class PointType>
std::vector<DGeo::Edge<PointType>> GridData3d<PointType>::GetCellEdgesArray(int cellNumb) const
{
    std::vector<DGeo::Edge<PointType>> edge(4);
    return edge;
}