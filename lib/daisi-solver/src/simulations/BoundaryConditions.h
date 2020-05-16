#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
}
class PropertyCondition
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    std::vector<std::vector<double>> conditionProperties;
    std::vector<std::vector<double>> zArray;
    std::vector<std::string>         propertiesNames;

  public:
    std::vector<std::string>          GetPropertiesSimpleNames() const;
    std::vector<std::vector<double>>& Get_zArray();
    std::vector<double>               GetSimpleProperties() const;
    void SetSimpleProperties(std::vector<double> in);
    double GetPotentialOffset() const;

    double GetPotentialAmplitude() const;
    double GetPotentialFrequency() const;
    double GetPotentialInitialPhase() const;

    double GetPotential(double t) const;
    double GetPotential(double t, double z, int& status) const;

    double GetPotentialOffset(double z, int& status) const;

    double GetPotentialAmplitude(double z, int& status) const;

    double GetPotentialFrequency(double z, int& status) const;

    double GetPotentialInitialPhase(double z, int& status) const;

    int    Typeflag;
    int    manualLineTag;
    int    AttachedElectrodeNumber;
    double manualLinePos;
    void SetPropertiesFromFile(std::string file);
    std::string      type;
    std::vector<int> boundariesList;
    PropertyCondition(std::string typeIn, int boundaryTypeFlag);
    PropertyCondition();
};
class BoundaryConditions
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  private:
    friend class boost::serialization::access;
    std::vector<int>               DefaultConditionsList;
    std::vector<PropertyCondition> PropertyConditionList;

  public:
    int         PropertyConditionListSize() const;
    std::string GetConditionPropertyType(int i) const;
    std::vector<double> GetConditionPropertiesSimple(int i) const;
    std::vector<std::string> GetConditionPropertiesSimpleNames(int i) const;

    void   clear();
    double GetPotential(int i, double t) const;
    double GetPotential(int i, double t, double z, int& status) const;
    double GetPotentialOffset(int i) const;
    double GetPotentialAmplitude(int i) const;
    int GetElectrodeNumber(int i) const;
    double GetPotentialFrequency(int i) const;
    double GetPotentialInitialPhase(int i) const;
    double GetPotentialOffset(int i, double z, int& status) const;
    double GetPotentialAmplitude(int i, double z, int& status) const;
    double GetPotentialFrequency(int i, double z, int& status) const;
    double GetPotentialInitialPhase(int i, double z, int& status) const;
    std::vector<double> GetPropertyConditionManualRestictions(int i);
    void SetPropertyConditionManualRestictions(int i, std::vector<double> params);
    void SetConditionProperties(int i, std::vector<double> cond);
    void SetDefaultConditionsList(const std::vector<int>& in);
    std::vector<int>               GetDefaultConditionsList() const;
    std::vector<PropertyCondition> GetPropertyConditionsList() const;
    void AddDefaultConditionsList(int i);
    void SetPropertyConditionsBoundariesList(int i, const std::vector<int>& in);
    void AddPropertyCondition(std::string type, int boundaryTypeFlag);
    std::vector<int> GetPropertyConditionsBoundariesList(int i) const;
    int GetPropertyConditionTypeFlag(int i) const;
    int  GetNumberProperties() const;
    void SetConditionPropertiesFromFile(int i, std::string cond);
    template <class PointType>
    bool CheckBoundaryConditionIntersection(const DGeo::Edge<PointType>& edge, int i);
};

#endif