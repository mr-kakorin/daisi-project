#include "GridData.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "ParticleShape2d.h"

template class GridData2d<double>;
template class GridData2d<float>;

template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_Ex()
{
    return this->E[0];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_Ey()
{
    return this->E[1];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_ExCol()
{
    return this->ECol[0];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_EyCol()
{
    return this->ECol[1];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_ExA()
{
    return this->EA[0];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Get_EyA()
{
    return this->EA[1];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Getx()
{
    return this->X[0];
};
template <class PointType>
std::vector<PointType>& GridData2d<PointType>::Gety()
{
    return this->X[1];
};

template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_Ex() const
{
    return this->E[0];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_Ey() const
{
    return this->E[1];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::GetECol() const
{
    return this->ECol[0];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_EyCol() const
{
    return this->ECol[1];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_ExCol() const
{
    return this->ECol[0];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_ExA() const
{
    return this->EA[0];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Get_EyA() const
{
    return this->EA[1];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Getx() const
{
    return this->X[0];
};
template <class PointType>
const std::vector<PointType>& GridData2d<PointType>::Gety() const
{
    return this->X[1];
};

template <class PointType>
void GridData2d<PointType>::GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                                    int PlotTypeFlag) const
{
    if (flag == flagStringsSolver::PlotFlags2d[4])
        Array[0] = (void*)(&this->rho[0][0]);

    size        = int(this->E[0].size());
    sizeElement = sizeof(this->E[0][0]);

    if (PlotTypeFlag == 0)
    {
        if (flag == flagStringsSolver::PlotFlags2d[0])
            Array[0] = (void*)(&this->E[0][0]);

        if (flag == flagStringsSolver::PlotFlags2d[1])
            Array[0] = (void*)(&this->E[1][0]);

        if (flag == flagStringsSolver::PlotFlags2d[3])
            Array[0] = (void*)(&this->V[0]);
    }

    if (PlotTypeFlag == 1)
    {
        if (flag == flagStringsSolver::PlotFlags2d[0])
            Array[0] = (void*)(&this->ECol[0][0]);

        if (flag == flagStringsSolver::PlotFlags2d[1])
            Array[0] = (void*)(&this->ECol[1][0]);

        if (flag == flagStringsSolver::PlotFlags2d[3])
            Array[0] = (void*)(&this->VCharge[0]);
    }
};

template <class PointType>
float GridData2d<PointType>::interpolatePoint(double x1, double x2, double, std::string value,
                                              int PlotTypeFlag) const
{
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
    ParticleShapeCIC2dStatic<PointType>::Wcalculate(PointType(x1), PointType(x2), this->X[0],
                                                    this->X[1], W, i, this->CICArray[i]);

    PointType result = 0;

    if (value == flagStringsSolver::PlotFlags2d[4])
        ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->rho[0], result, i,
                                                              this->CICArray[i]);

    if (PlotTypeFlag == 0)
    {

        if (value == flagStringsSolver::PlotFlags2d[0])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[1])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[2])
        {
            PointType tmp1;
            PointType tmp2;

            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], tmp1, i,
                                                                  this->CICArray[i]);
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], tmp2, i,
                                                                  this->CICArray[i]);

            result = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
        }

        if (value == flagStringsSolver::PlotFlags2d[3])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->V, result, i,
                                                                  this->CICArray[i]);
    }

    if (PlotTypeFlag == 1)
    {

        if (value == flagStringsSolver::PlotFlags2d[0])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[0], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[1])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[1], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[2])
        {
            PointType tmp1;
            PointType tmp2;

            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[0], tmp1, i,
                                                                  this->CICArray[i]);
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[1], tmp2, i,
                                                                  this->CICArray[i]);

            result = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
        }

        if (value == flagStringsSolver::PlotFlags2d[3])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->VCharge, result, i,
                                                                  this->CICArray[i]);
    }
    if (PlotTypeFlag == 2)
    {
        PointType result1 = 0;

        if (value == flagStringsSolver::PlotFlags2d[0])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[1])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[3])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->V, result, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[0])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[0], result1, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[1])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[1], result1, i,
                                                                  this->CICArray[i]);

        if (value == flagStringsSolver::PlotFlags2d[3])
            ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->VCharge, result1, i,
                                                                  this->CICArray[i]);

        result = result1 + result;
    }

    return result;
};

template <class PointType>
void GridData2d<PointType>::interpolatePoint(double x1, double x2, double x3, double& Exin,
                                             double& Eyin) const
{

    int i;
    for (i = 0; i < this->CICArray.size() - 1; i++)
    {
        if (ParticleShapeCIC2dStatic<PointType>::InCellWithEps(
                PointType(x1), PointType(x2), this->X[0], this->X[1], i, this->CICArray[i]))
            break;
    }
    if (i == this->CICArray.size() - 1)
    {
        Exin = 0;
        Eyin = 0;
        return;
    }

    std::vector<PointType> W(4);
    ParticleShapeCIC2dStatic<PointType>::Wcalculate(PointType(x1), PointType(x2), this->X[0],
                                                    this->X[1], W, i, this->CICArray[i]);

    PointType result  = 0;
    PointType result1 = 0;

    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[0], result, i,
                                                          this->CICArray[i]);
    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[0], result1, i,
                                                          this->CICArray[i]);

    Exin = result + result1;
    // Exin = result;
    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->E[1], result, i,
                                                          this->CICArray[i]);
    ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(W, this->ECol[1], result1, i,
                                                          this->CICArray[i]);

    // Eyin = result;
    Eyin = result + result1;
};

template <class PointType>
float GridData2d<PointType>::GetMaxSixe() const
{
    float Hmax = 0;
    float h;
    for (int i = 1; i < this->X[0].size(); i++)
    {

        h = std::abs(this->X[1][i] - this->X[1][i - 1]);

        if (h > Hmax)
            Hmax = h;
    }
    return std::max(float(std::abs(this->X[0][1] - this->X[0][0])), float(Hmax));
};

template <class PointType>
int GridData2d<PointType>::InCell(double x1, double x2, double x3) const
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
};

template <class PointType>
std::vector<DGeo::Edge<PointType>> GridData2d<PointType>::GetCellEdgesArray(int cellNumb) const
{
    std::vector<DGeo::Edge<PointType>> edge(4);

    DGeo::Point<PointType> p[4];

    p[0].x = this->X[0][cellNumb];
    p[0].y = this->X[1][cellNumb];
    p[0].z = 0;

    p[1].x = this->X[0][cellNumb + 1];
    p[1].y = this->X[1][cellNumb];
    p[1].z = 0;

    p[2].x = this->X[0][cellNumb + 1];
    p[2].y = this->X[1][this->CICArray[cellNumb]];
    p[2].z = 0;

    p[3].x = this->X[0][cellNumb];
    p[3].y = this->X[1][this->CICArray[cellNumb]];
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
};