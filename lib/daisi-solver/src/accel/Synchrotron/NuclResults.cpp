#include "Results.h"
#include "SynchrotronDevice.h"

void SimulationDataAccel::SetDataNucl(arma::mat& dx, arma::mat& dy, arma::mat& dx0, arma::mat& dy0, arma::vec& dIx,
                                      arma::vec& dIy, arma::vec& sx, arma::vec& sy)
{
    YData.resize(8);
    YData[0].resize(dx.n_rows);
    for (int i = 0; i < dx.n_rows; i++)
    {
        YData[0][i].resize(dx.n_cols);
        for (int j         = 0; j < dx.n_cols; j++)
            YData[0][i][j] = dx(i, j);
    };

    YData[1].resize(dy.n_rows);
    for (int i = 0; i < dy.n_rows; i++)
    {
        YData[1][i].resize(dy.n_cols);
        for (int j         = 0; j < dy.n_cols; j++)
            YData[1][i][j] = dy(i, j);
    };

    YData[2].resize(dx0.n_rows);
    for (int i = 0; i < dx0.n_rows; i++)
    {
        YData[2][i].resize(dx0.n_cols);
        for (int j         = 0; j < dx0.n_cols; j++)
            YData[2][i][j] = dx0(i, j);
    };

    YData[3].resize(dy0.n_rows);
    for (int i = 0; i < dy0.n_rows; i++)
    {
        YData[3][i].resize(dy0.n_cols);
        for (int j         = 0; j < dy0.n_cols; j++)
            YData[3][i][j] = dy0(i, j);
    };

    YData[4].resize(1);
    YData[4][0].resize(dIx.n_rows);
    for (int i         = 0; i < dIx.n_rows; i++)
        YData[4][0][i] = dIx(i);

    YData[5].resize(1);
    YData[5][0].resize(dIy.n_rows);
    for (int i         = 0; i < dIy.n_rows; i++)
        YData[5][0][i] = dIy(i);

    YData[6].resize(1);
    YData[6][0].resize(sx.n_rows);
    for (int i         = 0; i < sx.n_rows; i++)
        YData[6][0][i] = sx(i);

    YData[7].resize(1);
    YData[7][0].resize(sy.n_rows);
    for (int i         = 0; i < sy.n_rows; i++)
        YData[7][0][i] = sy(i);
};

void SimulationDataAccel::addAccelElemetsDescription(std::shared_ptr<SynchrotronDevice>& device)
{
    device->GetAccelElemetsDescription(props, names);
}
void SimulationDataAccel::initSynchrotron(int dataSize, int Ncircle, const std::vector<int>& ind,
                                          const std::vector<int>& indP)
{
    madXData.resize(Ncircle);
    for (int i = 0; i < Ncircle; i++)
    {
        madXData[i].resize(ind.size());
        for (int j = 0; j < ind.size(); j++)
        {
            madXData[i][j].resize(indP.size());
            for (int k = 0; k < indP.size(); k++)
                madXData[i][j][k].resize(6);
        }
    }
};
void SimulationDataAccel::setDataSynchrotron(int i, int circle, const std::vector<int>& ind,
                                             const std::vector<int>& indP, std::vector<void*> dataIn)
{
    int sizeElement = 8;
    int j;
    for (j = 0; j < ind.size(); j++)
    {
        if (i == ind[j])
            break;
    };
    if (j == ind.size())
        return;

    double tmp;

    for (int k = 0; k < indP.size(); k++)
    {
        for (int k1 = 0; k1 < 6; k1++)
        {
            tmp                        = *((double*)((char*)dataIn[k1] + indP[k] * sizeElement));
            madXData[circle][j][k][k1] = tmp;
        }
    }
};
