#include "ParticleSource.h"
#include "BoundaryContainer2d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "Geom.h"
#include "GridData.h"
#include "geomTools.h"
#include <Constants.h>
#include <sstream>

template <class PointType>
void ParticleSource2d<PointType>::SetFlowCurrent(double res)
{
    sourceCurrent = res;
}

template <class PointType>
bool ParticleSource2d<PointType>::EdgeCmp(DGeo::Edge<PointType>& Edge1,
                                          DGeo::Edge<PointType>& Edge2)
{
    double r1 = Edge1.Middle().Radius();
    double r2 = Edge2.Middle().Radius();
    return r1 < r2;
}

template <class PointType>
void ParticleSource2d<PointType>::setErAverage(double in)
{
    ErAverage = in;
}

template <class PointType>
double ParticleSource2d<PointType>::getErAverage()
{
    return ErAverage;
}

template <class PointType>
std::vector<DGeo::Point<double>> ParticleSource2d<PointType>::GetSufacePoints()
{
    std::vector<DGeo::Point<double>> x;
    for (int i = 0; i < sourceSurface.size(); i++)
    {
        DGeo::Point<PointType> tmp = sourceSurface[i].extractingEdge->Middle();
        DGeo::Point<double>    tmpx;
        tmpx.x = tmp.x;
        tmpx.y = tmp.y;
        tmpx.z = tmp.z;
        x.push_back(tmpx);
    }
    return x;
}

template <class PointType>
double ParticleSource2d<PointType>::length()
{
    return sourceSurface.back().curveLength +
           sourceSurface.back().extractingEdge->length() / 2; // for 2daxs case
}

template <class PointType>
template <class Archive>
void ParticleSource2d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& sourceSurface;
    ar& polinom;
}

template <class PointType>
template <class Archive>
void ParticleSource2d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& sourceSurface;
    ar& polinom;
    sourceCurrent = 0;
}

template<typename T>
void write_out( const std::vector<T> &to_out, const std::string &out_name ) {
	ofstream out_stream;
	out_stream.open( out_name );
	for (int i = 0; i < to_out.size(); ++i)
		out_stream << std::setprecision(16) << to_out[i] << std::endl;
	out_stream.close();
}

template <class PointType>
std::vector<PointType> ParticleSource2d<PointType>::GetParticle(PointType L1, PointType L2,
                                                                int flag)
{
    std::vector<PointType> out(7);
//    std::vector<double> points11, points22, points33, points44;
//    points11.reserve(sourceSurface.size());
//	points22.reserve(sourceSurface.size());
//	points33.reserve(sourceSurface.size());
//	points44.reserve(sourceSurface.size());
//    for(int i=0; i< sourceSurface.size(); ++i){
//	    points11.push_back(sourceSurface[i].extractingEdge->point1.x);
//	    points22.push_back(sourceSurface[i].extractingEdge->point1.y);
//	    points33.push_back(sourceSurface[i].extractingEdge->point2.x);
//	    points44.push_back(sourceSurface[i].extractingEdge->point2.y);
//    }
//
//	write_out<double>(points11, "./11.csv");
//	write_out<double>(points22, "./22.csv");
//	write_out<double>(points33, "./33.csv");
//	write_out<double>(points44, "./44.csv");
    PointType              lMiddle = (L1 + L2) / 2;
    int                    j       = 0;

    int j1 = sourceSurface.size() - 1;

    int middle;

    while (1)
    {
        if (j + 1 == j1 || j == j1)
            break;

        middle = (j + j1) / 2;

        if (lMiddle >= sourceSurface[middle].curveLength)
            j = middle;

        if (lMiddle <= sourceSurface[middle].curveLength)
            j1 = middle;
    }
    /*for (int i = lastIndex; i < sourceSurface.size(); i++)
    {
            if (sourceSurface[i - 1].curveLength <= lMiddle && sourceSurface[i].curveLength >=
    lMiddle)
            {
                    j = i - 1;
                    if (std::abs(sourceSurface[i - 1].curveLength - lMiddle) >
    std::abs(sourceSurface[i].curveLength - lMiddle))
                            j = i;
                    break;
            };
    };
    lastIndex = j;
    if (lastIndex == 0)
            lastIndex = 1;*/

    out[0]    = sourceSurface[j].extractingEdge->Middle().x;
    out[1]    = sourceSurface[j].extractingEdge->Middle().y;
    double jC = std::min(double(sourceSurface[j].currentDensity),
                         double(sourceSurface[j].maximalCurrentDensity));

    if (flag == 0)
        out[2] = (L2 - L1) * jC; // for 2d case
    if (flag == 1)
        out[2] = (L2 - L1) * jC * 2 * PI() * sourceSurface[j].extractingEdge->Middle().x;

    PointType tmp;
    PointType tmp1;

    Dmath::Cartesian2Polar(PointType(sourceSurface[j].normalX), PointType(sourceSurface[j].normalY),
                           tmp, tmp1);

    out[3] = PI() - tmp1;
    //	out[3] = sourceSurface[j].extractingEdge->alpha();

    out[4] = sourceSurface[j].alphaNormal;
	auto ff =
    // out[4] = out[0] + sourceSurface[j].normalX;
    out[5] = out[1] + sourceSurface[j].normalY;
    out[6] = PointType(sourceSurface[j].cellNumber);

    return out;
}

template <class PointType>
bool ParticleSource2d<PointType>::GetParticleOptimized(PointType L1, PointType L2, int flag, PointType*const& out, double const timestep)
{
	PointType              lMiddle = (L1 + L2) / 2;
	int                    j       = 0;

	int j1 = sourceSurface.size() - 1;

	int middle;

	while (true)
	{
		if (j + 1 == j1 || j == j1)
			break;

		middle = (j + j1) / 2;

		if (lMiddle >= sourceSurface[middle].curveLength)
			j = middle;

		if (lMiddle <= sourceSurface[middle].curveLength)
			j1 = middle;
	};

	int start = -1;
	int end = -1;
	for ( int i = 0; i < sourceSurface.size(); ++i)
	{
		if ( start == -1 && sourceSurface[i].curveLength >= L1 )
			start = i;
		if ( end == -1 && sourceSurface[i].curveLength >= L2 )
			end = i;
		if ( start != -1 && end != -1 )
			break;
	}

	out[2] = 0;
	for ( int i = start; i < end; ++i )
	{
		double jC = std::min(double(sourceSurface[i].currentDensity),
		                     double(sourceSurface[i].maximalCurrentDensity));

		if (flag == 0)
			out[2] += sourceSurface[i].extractingEdge->length() * jC; // for 2d case
		if (flag == 1)
			out[2] += sourceSurface[i].extractingEdge->length() * jC * 2 * PI() * sourceSurface[i].extractingEdge->Middle().x;
	}
	out[2] += sourceSurface[j].accumulatedCurrentDensity;
	std::stringstream ss;
	ss << j << "," << out[2] <<"," << timestep <<","<< timestep * out[2] <<","<< std::abs( ELECTRON_CHARGE() ) << "," <<std::endl;
	std::cout << ss.str();
	if ( std::isnan( out[2] ) ) {
		double p = sourceSurface[j].accumulatedCurrentDensity;
		double ppp = out[2];
	}


	if ( timestep * out[2] < std::abs( ELECTRON_CHARGE()) )
	{
		sourceSurface[j].accumulatedCurrentDensity = out[2];
		out[2] = 0;
	}
	else
	{
		sourceSurface[j].accumulatedCurrentDensity = 0;
	}
//	std::stringstream log;
//	log << out[2] << " vs " << ( (L2 - L1) * std::min(double(sourceSurface[j].currentDensity),double(sourceSurface[j].maximalCurrentDensity)) * 2 * PI() * sourceSurface[j].extractingEdge->Middle().x ) << " "
//	<< out[2] / ( (L2 - L1) * std::min(double(sourceSurface[j].currentDensity),double(sourceSurface[j].maximalCurrentDensity)) * 2 * PI() * sourceSurface[j].extractingEdge->Middle().x ) * 100 <<  "% new from old."  << std::endl;
//	std::cout << log.str();

	out[0]    = sourceSurface[j].extractingEdge->Middle().x;
	out[1]    = sourceSurface[j].extractingEdge->Middle().y;
//	double jC = std::min(double(sourceSurface[j].currentDensity),
//	                     double(sourceSurface[j].maximalCurrentDensity));
//
//	if (flag == 0)
//		out[2] = (L2 - L1) * jC; // for 2d case
//	if (flag == 1)
//		out[2] = (L2 - L1) * jC * 2 * PI() * sourceSurface[j].extractingEdge->Middle().x;

	PointType tmp;
	PointType tmp1;

	Dmath::Cartesian2Polar(PointType(sourceSurface[j].normalX), PointType(sourceSurface[j].normalY),
	                       tmp, tmp1);

	out[3] = PI() - tmp1;
	//	out[3] = sourceSurface[j].extractingEdge->alpha();

	out[4] = sourceSurface[j].alphaNormal;
	// out[4] = out[0] + sourceSurface[j].normalX;
	out[5] = out[1] + sourceSurface[j].normalY;
	out[6] = PointType(sourceSurface[j].cellNumber);

	auto x1 = sourceSurface[j].extractingEdge->point1.x;
	auto y1 = sourceSurface[j].extractingEdge->point1.y;

	auto x2 = sourceSurface[j].extractingEdge->point2.x;
	auto y2 = sourceSurface[j].extractingEdge->point2.y;
	PointType len = sqrt((y2 - y1)*(y2 - y1)+(x2 - x1)*(x2 - x1));
	out[7] = (y2-y1) / len;
	out[8] = -(x2-x1) / len;
	out[9] = 0;
	return out;
}

template <class PointType>
void ParticleSource2d<PointType>::resetSearch()
{
    lastIndex = 1;
}

template <class PointType>
double ParticleSource2d<PointType>::GetEmissionCurrent(int flag)
{
    if (sourceCurrent != 0 && !std::isnan(sourceCurrent) && !std::isinf(sourceCurrent))
        return sourceCurrent;

    double current = 0;
    if (flag == 0) // for 2d case
    {
        for (int i = 0; i < sourceSurface.size(); i++)
        {
            double jj =
                std::min(sourceSurface[i].currentDensity, sourceSurface[i].maximalCurrentDensity);
            current = current + sourceSurface[i].extractingEdge->length() * jj;
        }
    }
    if (flag == 1) // for 2daxs case
    {
        for (int i = 0; i < sourceSurface.size(); i++)
        {
            double jj =
                std::min(sourceSurface[i].currentDensity, sourceSurface[i].maximalCurrentDensity);
            current = current +
                      sourceSurface[i].extractingEdge->length() * jj * 2 * PI() *
                          sourceSurface[i].extractingEdge->Middle().x;
        }
    }
    return current;
}

template <class PointType>
template <class gridDataType>
void ParticleSource2d<PointType>::InitEmissionBoundary(
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<gridDataType>& grid, std::vector<double> parametersIn,
    std::string& errorMsg)
{
    sourceCurrent = 0;
    sourceSurface.clear();
    BoundaryContainer2d<PointType> boundary     = *boundaryIn[0];
    int                            flagBoundary = 1;

    std::vector<int> list;
    for (int i = 0; i < boundaryIn.size(); i++)
        list.push_back(i);

    std::vector<DGeo::Edge<PointType>> EdgesDataTmp;
    mergeSortResize(parametersIn[0], list, boundaryIn, EdgesDataTmp, errorMsg);

    if (errorMsg.size())
        return;

    DGeo::Point<PointType> pTest;
    pTest.x = parametersIn[1];
    pTest.y = parametersIn[2];
    pTest.z = parametersIn[3];

    int    nearestEdge = 0;
    double minDist     = pTest.Dist2Point(EdgesDataTmp[0].point1);
    for (int i = 1; i < EdgesDataTmp.size(); i++)
    {
        double dist = pTest.Dist2Point(EdgesDataTmp[i].point1);
        if (dist < minDist)
            nearestEdge = i;
    }

    float HP = 2 * grid->GetMaxSixe();

    DGeo::Point<PointType> ptmp1;
    DGeo::Point<PointType> ptmp2;

    ptmp1 = EdgesDataTmp[0].GetNormalPoint1(HP);
    ptmp1 = EdgesDataTmp.back().GetNormalPoint1(HP);

    ptmp1 = EdgesDataTmp[nearestEdge].GetNormalPoint1(HP);
    ptmp2 = EdgesDataTmp[nearestEdge].GetNormalPoint2(HP);

    int flagNormal;
    if (pTest.Dist2Point(ptmp1) < pTest.Dist2Point(ptmp2))
        flagNormal = 1;
    else
        flagNormal = 2;

    double                 curveLength = 0;
    DGeo::Point<PointType> ptmp;
    DGeo::Point<PointType> ptmpTest;

    for (int i = 0; i < EdgesDataTmp.size(); i++)
    {
        double alphaNormal;
        if (2 == flagNormal)
        {
            alphaNormal = EdgesDataTmp[i].alpha() + PI() / 2;
            //	ptmp = EdgesDataTmp[i].GetNormalPoint1(HP);
            // ptmpTest = EdgesDataTmp[i].GetNormalPoint2(HP);
        }
        else
        {
            alphaNormal = EdgesDataTmp[i].alpha() - PI() / 2;
            //	ptmp = EdgesDataTmp[i].GetNormalPoint2(HP);
            //	ptmpTest = EdgesDataTmp[i].GetNormalPoint1(HP);
        }

        double r     = sqrt(ptmp.x * ptmp.x + ptmp.y * ptmp.y);
        double rTest = sqrt(ptmpTest.x * ptmpTest.x + ptmpTest.y * ptmpTest.y);

        if (rTest > r)
        {
            int tt = 0;
        }

        int cellNumber = grid->InCell(EdgesDataTmp[i].Middle().x, EdgesDataTmp[i].Middle().y, 0);

        // sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesDataTmp[i], curveLength,
        // cellNumber, ptmp.x -
        // EdgesDataTmp[i].Middle().x, ptmp.y - EdgesDataTmp[i].Middle().y));
        sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesDataTmp[i], curveLength,
                                                              cellNumber, alphaNormal, flagNormal));

        curveLength = curveLength + EdgesDataTmp[i].length();
    }

    /*double length=0;
    for (int i = 0; i < size; i++)
    {
            length = length + EdgesData[i].length();
            //	currentDensityDistribution.push_back();
    };
    double avLength = length / normSize;
    double curveLength = 0;
    int cellNumber;
    DGeo::Point<PointType> ptmp;
    for (int i = 0; i < size; i++)
    {
            if (EdgesData[i].length() < avLength)
            {
                    //EdgesData[i]
                    ptmp = EdgesData[i].GetNormalPoint1(HP);

                    cellNumber = grid->InCell(EdgesData[i].Middle().x, EdgesData[i].Middle().y, 0);

                    if (grid->InCell(ptmp.x, ptmp.y, 0) == -1)
                            sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesData[i],
    curveLength, cellNumber,
    ptmp.x - EdgesData[i].Middle().x, ptmp.y - EdgesData[i].Middle().y)); else
                    {
                            ptmp = EdgesData[i].GetNormalPoint2(HP);
                            sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesData[i],
    curveLength, cellNumber,
    ptmp.x - EdgesData[i].Middle().x, ptmp.y - EdgesData[i].Middle().y));
                    }
                    curveLength = curveLength + EdgesData[i].length();
            }
            else
            {
                    int resize = ceil(EdgesData[i].length() / avLength);
                    std::vector<DGeo::Edge <PointType>> resArr = EdgesData[i].resize(resize);
                    for (int j = 0; j < resize; j++)
                    {
                            ptmp = resArr[j].GetNormalPoint1(HP);

                            cellNumber = grid->InCell(resArr[j].Middle().x, resArr[j].Middle().y,
    0);

                            if (grid->InCell(ptmp.x, ptmp.y, 0) == -1)
                                    sourceSurface.push_back(ParticleSourceEdge<PointType>(resArr[j],
    curveLength,
    cellNumber, ptmp.x - resArr[j].Middle().x, ptmp.y - resArr[j].Middle().y)); else
                            {
                                    ptmp = resArr[j].GetNormalPoint2(HP);
                                    sourceSurface.push_back(ParticleSourceEdge<PointType>(resArr[j],
    curveLength,
    cellNumber, ptmp.x - resArr[j].Middle().x, ptmp.y - resArr[j].Middle().y));
                            }
                            curveLength = curveLength + resArr[j].length();
                    };
            };
    };

    polinom.resize(11);

    for (int i = 0; i < polinom.size(); i++)
            polinom[i] = 0;*/
}

template <class PointType>
void ParticleSource2d<PointType>::InitEmissionBoundary(
    std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
    const std::shared_ptr<GridData2dpolar<PointType>>& grid, std::vector<double> parametersIn,
    std::string& errorMsg)
{
    sourceCurrent = 0;
    sourceSurface.clear();
    BoundaryContainer2d<PointType> boundary = *boundaryIn[0];

    int flagBoundary = 1;

    //	for (int i = 1; i < boundaryIn.size(); i++)
    //		flagBoundary = boundary.MergeWithSort(boundaryIn[i]);

    if (!flagBoundary)
    {
        int error = 1;
    };

    int   size     = boundary.EdgesData.size();
    int   normSize = 400;
    float HP       = 2 * grid->GetMaxSixe();

    std::vector<DGeo::Edge<PointType>> EdgesData = boundary.EdgesData;

    //	if (EdgesData.size()>1)
    //	std::sort(EdgesData.begin(), EdgesData.end(), EdgeCmp);

    // if (size >= normSize)
    double length = 0;
    for (int i = 0; i < size; i++)
    {
        length = length + EdgesData[i].length();
        //	currentDensityDistribution.push_back();
    };
    double                 avLength    = length / normSize;
    double                 curveLength = 0;
    int                    cellNumber;
    DGeo::Point<PointType> ptmp;

    PointType rTmp;
    PointType phiTmp;

    HP = HP / 10;
    for (int i = 0; i < size; i++)
    {
        if (EdgesData[i].length() < avLength)
        {
            // EdgesData[i]
            ptmp = EdgesData[i].GetNormalPoint1(HP);

            Dmath::Cartesian2Polar(EdgesData[i].Middle().x, EdgesData[i].Middle().y, rTmp, phiTmp);

            cellNumber = grid->InCellWithEps(rTmp, phiTmp, 0);

            Dmath::Cartesian2Polar(ptmp.x, ptmp.y, rTmp, phiTmp);

            if (grid->InCellWithEps(rTmp, phiTmp, 0) == -1)
                sourceSurface.push_back(ParticleSourceEdge<PointType>(
                    EdgesData[i], curveLength, cellNumber, ptmp.x - EdgesData[i].Middle().x,
                    ptmp.y - EdgesData[i].Middle().y));
            else
            {
                ptmp = EdgesData[i].GetNormalPoint2(HP);
                sourceSurface.push_back(ParticleSourceEdge<PointType>(
                    EdgesData[i], curveLength, cellNumber, ptmp.x - EdgesData[i].Middle().x,
                    ptmp.y - EdgesData[i].Middle().y));
            }
            curveLength = curveLength + EdgesData[i].length();
        }
        else
        {
            int                                resize = ceil(EdgesData[i].length() / avLength);
            std::vector<DGeo::Edge<PointType>> resArr = EdgesData[i].resize(resize);
            for (int j = 0; j < resize; j++)
            {
                int flag = 0;
                for (int pt = 1; pt < 100; pt++)
                {
                    ptmp = resArr[j].GetNormalPoint1(HP * 0.01 * pt);

                    Dmath::Cartesian2Polar(resArr[j].Middle().x, resArr[j].Middle().y, rTmp,
                                           phiTmp);

                    cellNumber = grid->InCellWithEps(rTmp, phiTmp, 0);

                    Dmath::Cartesian2Polar(ptmp.x, ptmp.y, rTmp, phiTmp);

                    if (grid->InCellWithEps(rTmp, phiTmp, 0) == -1)
                    {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0)
                {
                    ptmp = resArr[j].GetNormalPoint1(HP);
                    sourceSurface.push_back(ParticleSourceEdge<PointType>(
                        resArr[j], curveLength, cellNumber, ptmp.x - resArr[j].Middle().x,
                        ptmp.y - resArr[j].Middle().y));
                }
                else
                {
                    ptmp = resArr[j].GetNormalPoint2(HP);
                    sourceSurface.push_back(ParticleSourceEdge<PointType>(
                        resArr[j], curveLength, cellNumber, ptmp.x - resArr[j].Middle().x,
                        ptmp.y - resArr[j].Middle().y));
                }
                curveLength = curveLength + resArr[j].length();
            }
        }
    }

    polinom.resize(11);

    for (int i     = 0; i < polinom.size(); i++)
        polinom[i] = 0;
}
/*template<class PointType>
void ParticleSource2d<PointType> ::InitEmissionBoundary(std::vector<BoundaryContainer2d <PointType>>
boundaryIn,const std::shared_ptr<GridData2daxs<PointType>>& grid)
{
        sourceSurface.clear();
        BoundaryContainer2d <PointType> boundary = boundaryIn[0];

        if (boundary.EdgesData[0].Middle().y > boundary.EdgesData.back().Middle().y)
                boundary.Reverse();

        int flagBoundary = 1;

        for (int i = 1; i < boundaryIn.size(); i++)
                flagBoundary = boundary.MergeWithSort(boundaryIn[i]);

        if (!flagBoundary)
        {
                int error = 1;
        };

        int size = boundary.ContainerSize;
        int normSize = 15000;

        std::vector<DGeo::Edge<PointType>> EdgesData = boundary.EdgesData;

        //	if (EdgesData.size()>1)
        //	std::sort(EdgesData.begin(), EdgesData.end(), EdgeCmp);

        //if (size >= normSize)

        float HP = 2 * grid->GetMaxSixe();

        double length = 0;
        for (int i = 0; i < size; i++)
        {
                length = length + EdgesData[i].length();
                //	currentDensityDistribution.push_back();
        };
        double avLength = length / normSize;
        double curveLength = 0;
        int cellNumber;
        DGeo::Point<PointType> ptmp;
        for (int i = 0; i < size; i++)
        {
                if (EdgesData[i].length() < avLength)
                {
                        //EdgesData[i]

                        ptmp = EdgesData[i].GetNormalPoint1(HP);

                        cellNumber = grid->InCell(EdgesData[i].Middle().x, EdgesData[i].Middle().y,
0);

                        if (grid->InCell(ptmp.x, ptmp.y, 0) == -1)
                                sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesData[i],
curveLength,
cellNumber, ptmp.x - EdgesData[i].Middle().x, ptmp.y - EdgesData[i].Middle().y)); else
                        {
                                ptmp = EdgesData[i].GetNormalPoint2(HP);
                                sourceSurface.push_back(ParticleSourceEdge<PointType>(EdgesData[i],
curveLength,
cellNumber, ptmp.x - EdgesData[i].Middle().x, ptmp.y - EdgesData[i].Middle().y));
                        }

                        curveLength = curveLength + EdgesData[i].length();
                }
                else
                {
                        int resize = ceil(EdgesData[i].length() / avLength);
                        std::vector<DGeo::Edge <PointType>> resArr = EdgesData[i].resize(resize);
                        for (int j = 0; j < resize; j++)
                        {
                                ptmp = resArr[j].GetNormalPoint1(HP);

                                cellNumber = grid->InCell(resArr[j].Middle().x,
resArr[j].Middle().y, 0);

                                if (grid->InCell(ptmp.x, ptmp.y, 0) == -1)
                                        sourceSurface.push_back(ParticleSourceEdge<PointType>(resArr[j],
curveLength,
cellNumber, ptmp.x - resArr[j].Middle().x, ptmp.y - resArr[j].Middle().y)); else
                                {
                                        ptmp = resArr[j].GetNormalPoint2(HP);
                                        sourceSurface.push_back(ParticleSourceEdge<PointType>(resArr[j],
curveLength,
cellNumber, ptmp.x - resArr[j].Middle().x, ptmp.y - resArr[j].Middle().y));
                                }

                                curveLength = curveLength + resArr[j].length();
                        };
                };
        };

        polinom.resize(11);

        for (int i = 0; i < polinom.size(); i++)
                polinom[i] = 0;
}*/

template <class PointType>
void ParticleSource2d<PointType>::InitEmissionBoundary(
    const std::shared_ptr<ElectrodeCurrent<PointType>>& electrodeIn,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, std::vector<double> parametersIn,
    std::string& errorMsg)
{
    electrode = electrodeIn;
    BoundaryContainer2d<PointType> boundary;
    boundary.EdgesData = electrodeIn->ElectrodeEdges;
    InitEmissionBoundary({std::shared_ptr<BoundaryContainer2d<PointType>>(&boundary)}, gridData,
                         parametersIn, errorMsg);
}

template <class PointType>
std::vector<std::vector<float>> ParticleSource2d<PointType>::GetCurrentDensityDistribution()
{
    std::vector<float> v1;
    std::vector<float> v2;
    for (int i = 0; i < sourceSurface.size(); i++)
    {
        v1.push_back(sourceSurface[i].curveLength);
        v2.push_back(
            std::min(sourceSurface[i].currentDensity, sourceSurface[i].maximalCurrentDensity));
    }
    std::vector<std::vector<float>> result;
    result.push_back(v1);
    result.push_back(v2);
    return result;
}

template <class PointType>
std::vector<std::vector<double>> ParticleSource2d<PointType>::GetEmitterField()
{
    std::vector<double> v1;
    std::vector<double> v2;
    for (int i = 0; i < sourceSurface.size(); i++)
    {
        v1.push_back(sourceSurface[i].curveLength);
        v2.push_back(sourceSurface[i].E);
    }
    std::vector<std::vector<double>> result;
    result.push_back(v1);
    result.push_back(v2);
    return result;
}

template <class PointType>
void ParticleSource2d<PointType>::SetZCoordinate(double Z)
{
    for (int i = 0; i < sourceSurface.size(); i++)
    {
        sourceSurface[i].extractingEdge->point1.z = Z;
        sourceSurface[i].extractingEdge->point2.z = Z;
    }
}

template <class PointType>
template <class Archive>
void ParticleSourceEdge<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& curveLength;
    ar& currentDensity;
    ar& extractingEdge;
    ar& alphaNormal;
    ar& flagNormal;
    ar& cellNumber;
    ar& maximalCurrentDensity;
}

template <class PointType>
template <class Archive>
void ParticleSourceEdge<PointType>::load(Archive& ar, const unsigned int)
{
    ar& curveLength;
    ar& currentDensity;
    ar& extractingEdge;
    ar& alphaNormal;
    ar& flagNormal;
    ar& cellNumber;
    ar& maximalCurrentDensity;
    E = 0;
}

template <class PointType>
ParticleSourceEdge<PointType>::ParticleSourceEdge(DGeo::Edge<PointType>& Edge,
                                                  double curveLengthPrev, int cell, double X,
                                                  double Y)
{
    extractingEdge        = std::shared_ptr<DGeo::Edge<PointType>>(new DGeo::Edge<PointType>());
    *extractingEdge       = Edge;
    currentDensity        = 0;
    maximalCurrentDensity = 1;
    curveLength           = curveLengthPrev + Edge.length() / 2;
    cellNumber            = cell;
    alphaNormal           = X;
    flagNormal            = Y;

	accumulatedCurrentDensity = 0;
}

template <class PointType>
ParticleSourceEdge<PointType>::ParticleSourceEdge(){

}

template class ParticleSource2d<float>;
template class ParticleSource2d<double>;

template void ParticleSource2d<double>::InitEmissionBoundary<GridData2daxs<double>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<double>>> boundaryIn,
        const std::shared_ptr<GridData2daxs<double>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSource2d<float>::InitEmissionBoundary<GridData2daxs<float>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<float>>> boundaryIn,
        const std::shared_ptr<GridData2daxs<float>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSource2d<double>::InitEmissionBoundary<GridData2d<double>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<double>>> boundaryIn,
        const std::shared_ptr<GridData2d<double>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSource2d<float>::InitEmissionBoundary<GridData2d<float>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<float>>> boundaryIn,
        const std::shared_ptr<GridData2d<float>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSource2d<double>::InitEmissionBoundary<GridData3d<double>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<double>>> boundaryIn,
        const std::shared_ptr<GridData3d<double>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSource2d<float>::InitEmissionBoundary<GridData3d<float>>(
        std::vector<std::shared_ptr<BoundaryContainer2d<float>>> boundaryIn,
        const std::shared_ptr<GridData3d<float>>& grid, std::vector<double> parametersIn,
        std::string& errorMsg);

template void ParticleSourceEdge<float>::load<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ParticleSourceEdge<double>::load<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version);

template void ParticleSourceEdge<double>::load<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ParticleSourceEdge<float>::load<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version);

template void
ParticleSource2d<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);
template void ParticleSource2d<double>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void
ParticleSource2d<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void ParticleSource2d<float>::save<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version) const;