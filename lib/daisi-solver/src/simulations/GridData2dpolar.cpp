#include "GridData.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "ParticleShape2d.h"


template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_ErCol()
{
    return this->ECol[0];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_EphiCol()
{
    return this->ECol[1];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_Er()
{
    return this->E[0];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_Ephi()
{
    return this->E[1];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_ErA()
{
    return this->EA[0];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Get_EphiA()
{
    return this->EA[1];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Getr()
{
    return this->X[0];
}

template <class PointType>
std::vector<PointType>& GridData2dpolar<PointType>::Getphi()
{
    return this->X[1];
}

template <class PointType>
float GridData2dpolar<PointType>::GetMaxSixe() const
{
    float Hmax = 0;
    float h;
    for (int i = 1; i < this->X[0].size(); i++)
    {
        h = std::abs(this->X[0][i] - this->X[0][i - 1]);

        if (h > Hmax && this->X[1][i] == this->X[1][i - 1])
            Hmax = h;
    }
    return Hmax;
}

/*template void GridData2dpolar<double>::setF(GridData2daxs<double>& gr);
template void GridData2dpolar<double>::setF(GridData2daxs<float>& gr);

template<class PointType>
template <class PointType1>
void GridData2dpolar<PointType>::setF(GridData2daxs<PointType1>& gr)
{
        /*	densityReset();
        for (int i = 0; i < this->X[0].size(); i++)
        {
        for (int j = 0; j < gr.r.size()-1; j++)
        {
        if (this->X[0][i] >= gr.this->X[0][j] && this->X[0][i] <= gr.this->X[0][j+1])
        {
        double w1 = (gr.this->X[0][j + 1] - this->X[0][i]) / (gr.this->X[0][j + 1] -
gr.this->X[0][j]); double w2 = 1 - w1;

        PointType tmp = gr.Getthis->rho()[j] * w1 + gr.Getthis->rho()[j + 1] * w2;
        this->rho[0][i] = tmp;
        break;
        };
        };
        };
};*/

template <class PointType>
int GridData2dpolar<PointType>::InCell(double x1, double x2, double x3) const
{
    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCell(PointType(x1), PointType(x2), this->X[0],
                                                        this->X[1], i, this->CICArray[i]))
            break;
    }
    if (i == this->CICArray.size() - 1)
        return -1;
    return i;
}

template <class PointType>
float GridData2dpolar<PointType>::interpolatePoint(double x1, double x2, double, std::string value,
                                                   int PlotTypeFlag) const
{
    Dmath::Cartesian2Polar(x1, x2, x1, x2);

    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCellWithEps(
                PointType(x1), PointType(x2), this->X[0], this->X[1], i, this->CICArray[i]))
            break;
    }
    if (i == this->CICArray.size() - 1)
        return NAN;

    std::vector<PointType> W(4);
    ParticleShapeCIC2dStatic<PointType>::WcalculatePolar(PointType(x1), PointType(x2), this->X[0],
                                                         this->X[1], W, i, this->CICArray[i]);

    PointType result = 0;

    if (value == flagStringsSolver::PlotFlags2dpolar[0])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], result, i,
                                                              this->CICArray[i]);

    if (value == flagStringsSolver::PlotFlags2dpolar[1])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], result, i,
                                                              this->CICArray[i]);

    if (value == flagStringsSolver::PlotFlags2daxs[2])
    {
        PointType tmp1;
        PointType tmp2;

        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], tmp1, i,
                                                              this->CICArray[i]);
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], tmp2, i,
                                                              this->CICArray[i]);

        result = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
    }

    if (value == flagStringsSolver::PlotFlags2dpolar[3])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->V, result, i,
                                                              this->CICArray[i]);

    if (value == flagStringsSolver::PlotFlags2dpolar[4])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->rho[0], result, i,
                                                              this->CICArray[i]);

    return result;
}

template <class PointType>
int GridData2dpolar<PointType>::InCellWithEps(double x1, double x2, double x3) const
{
    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCellWithEps(
                PointType(x1), PointType(x2), this->X[0], this->X[1], i, this->CICArray[i]))
            break;
    }
    if (i == this->CICArray.size() - 1)
        return -1;
    return i;
}

template <class PointType>
std::vector<DGeo::Edge<PointType>> GridData2dpolar<PointType>::GetCellEdgesArray(int cellNumb) const
{
    std::vector<DGeo::Edge<PointType>> edge(4);

    DGeo::Point<PointType> p[4];

    PointType xtmp, ytmp;

    Dmath::Polar2Cartesian(this->X[0][cellNumb], this->X[1][cellNumb], xtmp, ytmp);

    p[0].x = xtmp;
    p[0].y = ytmp;
    p[0].z = 0;

    Dmath::Polar2Cartesian(this->X[0][cellNumb + 1], this->X[1][cellNumb], xtmp, ytmp);

    p[1].x = xtmp;
    p[1].y = ytmp;
    p[1].z = 0;

    Dmath::Polar2Cartesian(this->X[0][cellNumb + 1], this->X[1][this->CICArray[cellNumb]], xtmp,
                           ytmp);

    p[2].x = xtmp;
    p[2].y = ytmp;
    p[2].z = 0;

    Dmath::Polar2Cartesian(this->X[0][cellNumb], this->X[1][this->CICArray[cellNumb]], xtmp, ytmp);

    p[3].x = xtmp;
    p[3].y = ytmp;
    p[3].z = 0;

    edge[0].point1 = p[0];
    edge[0].point2 = p[1];

    edge[1].point1 = p[1];
    edge[1].point2 = p[2];

    edge[2].point1 = p[2];
    edge[2].point2 = p[3];

    edge[3].point1 = p[3];
    edge[3].point2 = p[0];

    return edge;
}

template <class PointType>
void GridData2dpolar<PointType>::interpolatePoint(double x1, double x2, double x3, double& Exin,
                                                  double& Eyin) const
{
    Dmath::Cartesian2Polar(x1, x2, x1, x2);

    double ErTmp, EphiTmp;

    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCellWithEps(
                PointType(x1), PointType(x2), this->X[0], this->X[1], i, this->CICArray[i]))
            break;
    }
    if (i == this->CICArray.size() - 1)
    {
        ErTmp   = 0;
        EphiTmp = 0;
        return;
    }

    std::vector<PointType> W(4);
    ParticleShapeCIC2dStatic<PointType>::Wcalculate(PointType(x1), PointType(x2), this->X[0],
                                                    this->X[1], W, i, this->CICArray[i]);

    PointType result = 0;

    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], result, i,
                                                          this->CICArray[i]);
    ErTmp = result;
    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], result, i,
                                                          this->CICArray[i]);
    EphiTmp = result;

    Exin = ErTmp * std::cos(2 * PI() - x2) - EphiTmp * std::sin(2 * PI() - x2);

    Eyin = ErTmp * std::sin(2 * PI() - x2) + EphiTmp * std::cos(2 * PI() - x2);
}

template class GridData2dpolar<double>;
template class GridData2dpolar<float>;