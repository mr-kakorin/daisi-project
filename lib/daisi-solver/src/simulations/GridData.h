#ifndef GRIDDATA_H
#define GRIDDATA_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
};
namespace Dmath
{
class imat;
};
template <class PointType>
class NearCathodeVolume
{

  public:
    std::vector<DGeo::Edge<PointType>> Edge;

    double fieldPointsX[4];
    double fieldPointsY[4];
    double Areas[4];
    double normalX[4];
    double normalY[4];
    double volume;

    int flagCurrentLimited;

    bool InCell(PointType x1, PointType x2) const;

    NearCathodeVolume(DGeo::Edge<PointType> Edge1In, DGeo::Edge<PointType> Edge2In, int flag,
                      int flagEm);
};

template <class PointType>
class GridDataBase
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  protected:
    std::vector<std::vector<PointType>> X;
    std::vector<std::vector<PointType>> E;
    std::vector<std::vector<PointType>> EA;
    std::vector<std::vector<PointType>> ECol;

    std::vector<std::vector<PointType>> B;

    std::vector<std::vector<PointType>> rho;
    std::vector<PointType>              V;
    std::vector<PointType>              VA;
    std::vector<PointType>              VCharge;
    std::vector<PointType>              flagOut;
    std::vector<int>                    CICArray;

  public:
    void Init();
    void ZeroingFields();
    void ApplyTimeDepending(PointType frequency, PointType phase, PointType time,
                            std::vector<int>& nonZeros);
    void ApplyTimeDepending(PointType frequency, PointType phase, PointType time);
    void ApplyTimeDepending(const std::vector<double>& globalPar, double time);
    void Clear();
    void Init(std::vector<DGeo::Point<PointType>> serialMeshData, int problemType);
    void SetCells(Dmath::imat& flagMatrix, Dmath::imat& templNumb, int nVertexX, int nVertexY,
                  int problemType);

    void densityReset();

    std::vector<PointType>& Getrho();

    std::vector<PointType>& GetV();

    std::vector<PointType>& GetVCharge();
    std::vector<PointType>& GetVA();

    void InitParallel(int numThreads);

    std::vector<PointType>& Getrho(int thread);

    void Summrho();

  public:
    void GetDataIntFlag(void* Array[1], int& size, int& sizeElement, int flag,
                        int PlotTypeFlag) const;
    void ZeroingFieldsBase();
    GridDataBase();
    std::vector<std::vector<PointType>> rhoParallel;
    bool SearchIntersectionWithEdge(int cellNumb, DGeo::Edge<PointType> EdgeIntersection,
                                    DGeo::Point<PointType>* startPoint) const;
    virtual std::vector<DGeo::Edge<PointType>> GetCellEdgesArray(int cellNumb) const = 0;
};

template <class PointType>
class GridData1d
{
    friend class boost::serialization::access;

  public:
    std::vector<PointType> x;
    std::vector<PointType> flagOut;
    std::vector<PointType> E[0];
};
template <class PointType>
class GridData2d : public GridDataBase<PointType>
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };

  public:
    std::vector<PointType>& Get_Ex();
    std::vector<PointType>& Get_Ey();
    std::vector<PointType>& Get_ExCol();
    std::vector<PointType>& Get_EyCol();
    std::vector<PointType>& Get_ExA();
    std::vector<PointType>& Get_EyA();
    std::vector<PointType>& Getx();
    std::vector<PointType>& Gety();

    const std::vector<PointType>& Get_Ex() const;
    const std::vector<PointType>& Get_Ey() const;
    const std::vector<PointType>& GetECol() const;
    const std::vector<PointType>& Get_EyCol() const;
    const std::vector<PointType>& Get_ExCol() const;
    const std::vector<PointType>& Get_ExA() const;
    const std::vector<PointType>& Get_EyA() const;
    const std::vector<PointType>& Getx() const;
    const std::vector<PointType>& Gety() const;

    std::vector<DGeo::Edge<PointType>> GetCellEdgesArray(int cellNumb) const;

    int GetType() const
    {
        return 1;
    };
    PointType* GetX1Pointer() const
    {
        return (PointType*)(&this->X[0][0]);
    };
    PointType* GetX2Pointer() const
    {
        return (PointType*)(&this->X[1][0]);
    };
    float GetMaxSixe() const;

    std::vector<PointType> GetB() const
    {
        std::vector<PointType> r;
        return r;
    };
    void SetB(const std::vector<PointType>& BIn){};
    void GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                 int PlotTypeFlag) const;
    float interpolatePoint(double x1, double x2, double x3, std::string value,
                           int PlotTypeFlag) const;
    void interpolatePoint(double x1, double x2, double x3, double& Ex, double& Ey) const;
    int InCell(double x1, double x2, double x3) const;
};

template <class PointType>
class GridData3d : public GridDataBase<PointType>
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };

  public:
    std::vector<PointType>& Getx();
    std::vector<PointType>& Gety();
    std::vector<PointType>& Getz();
    std::vector<PointType>& Get_Ex();
    std::vector<PointType>& Get_Ey();
    std::vector<PointType>& Get_Ez();
    std::vector<PointType>& Get_ExA();
    std::vector<PointType>& Get_EyA();
    std::vector<PointType>& Get_EzA();
    std::vector<PointType>& Get_ExCol();
    std::vector<PointType>& Get_EyCol();
    std::vector<PointType>& Get_EzCol();
    std::vector<PointType>& GetBx();
    std::vector<PointType>& GetBy();
    std::vector<PointType>& GetBz();

    std::vector<int> CICArrayZ;

    std::vector<DGeo::Edge<PointType>> GetCellEdgesArray(int cellNumb) const;

    int GetType() const
    {
        return 4;
    };
    PointType* GetX1Pointer() const
    {
        return (PointType*)(&this->X[0][0]);
    };
    PointType* GetX2Pointer() const
    {
        return (PointType*)(&this->X[1][0]);
    };

    PointType* GetX3Pointer() const
    {
        return (PointType*)(&this->X[2][0]);
    };

    float GetMaxSixe() const;

    std::vector<PointType> GetB() const
    {
        std::vector<PointType> r;
        return r;
    };
    void SetB(const std::vector<PointType>& BIn){};
    void GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                 int PlotTypeFlag) const;
    float interpolatePoint(double x1, double x2, double x3, std::string value,
                           int PlotTypeFlag) const;
    void interpolatePoint(double x1, double x2, double x3, double& Ex, double& Ey) const;
    int InCell(double x1, double x2, double x3) const;
};

template <class PointType>
class GridData2daxs : public GridDataBase<PointType>
{
    friend class boost::serialization::access;

  public:
    std::vector<PointType>& Get_Er();
    std::vector<PointType>& Get_Ez();
    std::vector<PointType>& Get_ErCol();
    std::vector<PointType>& Get_EzCol();
    std::vector<PointType>& Get_ErA();
    std::vector<PointType>& Get_EzA();
    std::vector<PointType>& Getr();
    std::vector<PointType>& Getz();
    std::vector<PointType>& GetBphi();
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    int GetType() const
    {
        return 2;
    };
    std::vector<DGeo::Edge<PointType>> GetCellEdgesArray(int cellNumb) const;
    float GetMaxSixe() const;
    int InCell(double x1, double x2, double x3) const;

    PointType* GetX1Pointer() const
    {
        return (PointType*)(&this->X[0][0]);
    };
    PointType* GetX2Pointer() const
    {
        return (PointType*)(&this->X[1][0]);
    };
    void GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                 int PlotTypeFlag) const;
    double interpolatePoint(double x1, double x2, double, std::string value, int PlotTypeFlag) const;
    void interpolatePoint(double x1, double x2, double x3, double& Ex, double& Ey) const;
};
template <class PointType>
class GridData2dpolar : public GridDataBase<PointType>
{
    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<GridDataBase<PointType>>(*this);
    };

    std::vector<PointType>& Get_ErCol();
    std::vector<PointType>& Get_EphiCol();
    std::vector<PointType>& Get_Er();
    std::vector<PointType>& Get_Ephi();
    std::vector<PointType>& Get_ErA();
    std::vector<PointType>& Get_EphiA();
    std::vector<PointType>& Getr();
    std::vector<PointType>& Getphi();

    std::vector<DGeo::Edge<PointType>> GetCellEdgesArray(int cellNumb) const;
    int GetType() const
    {
        return 3;
    };
    PointType* GetX1Pointer() const
    {
        return (PointType*)(&this->X[0][0]);
    };
    PointType* GetX2Pointer() const
    {
        return (PointType*)(&this->X[1][0]);
    };

    float GetMaxSixe() const;
    int InCell(double x1, double x2, double x3) const;
    void GetData(void* Array[1], int& size, int& sizeElement, std::string flag,
                 int PlotTypeFlag) const {

    };

    std::vector<PointType> GetB()
    {
        std::vector<PointType> tmp;
        return tmp;
    };
    void SetB(const std::vector<PointType>& BIn){};

    float interpolatePoint(double x1, double x2, double, std::string value, int PlotTypeFlag) const;
    int InCellWithEps(double x1, double x2, double x3) const;
    void interpolatePoint(double x1, double x2, double x3, double& Ex, double& Ey) const;
};

#endif