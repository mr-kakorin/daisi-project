#include "BoundaryConditions.h"
#include "Geom.h"
#include "VTKIncludeSolver.h"
#include <common_tools/constants.h>

const static std::vector<std::string> flowBoundaryTypeNames = {
    "Transparency", "Reflection", "Backscattering", "Absolut Absorbtion"};
std::vector<std::vector<std::string>> properiesNames = {
    {"Transparency value"}, {}, {"alpha fraction", "beta fraction"}, {},
};
const static std::vector<std::string> fieldParameters = {
    std::string("potential offset, V"), std::string("potential amplitude, V"),
    std::string("frequency, Hz"), std::string("initial phase, rad")};

void BoundaryConditions::clear()
{
    DefaultConditionsList.clear();
    PropertyConditionList.clear();
}

double BoundaryConditions::GetPotential(int i, double t) const
{
    return PropertyConditionList[i].GetPotential(t);
};

double BoundaryConditions::GetPotential(int i, double t, double z, int& status) const
{
    return PropertyConditionList[i].GetPotential(t, z, status);
};
double BoundaryConditions::GetPotentialOffset(int i) const
{
    return PropertyConditionList[i].GetPotentialOffset();
};
double BoundaryConditions::GetPotentialAmplitude(int i) const
{
    return PropertyConditionList[i].GetPotentialAmplitude();
};

int BoundaryConditions::GetElectrodeNumber(int i) const
{
    // return PropertyConditionList[i].AttachedElectrodeNumber;
    return -1;
};

double BoundaryConditions::GetPotentialFrequency(int i) const
{
    return PropertyConditionList[i].GetPotentialFrequency();
};
double BoundaryConditions::GetPotentialInitialPhase(int i) const
{
    return PropertyConditionList[i].GetPotentialInitialPhase();
};

double BoundaryConditions::GetPotentialOffset(int i, double z, int& status) const
{
    return PropertyConditionList[i].GetPotentialOffset(z, status);
};
double BoundaryConditions::GetPotentialAmplitude(int i, double z, int& status) const
{
    return PropertyConditionList[i].GetPotentialAmplitude(z, status);
};

double BoundaryConditions::GetPotentialFrequency(int i, double z, int& status) const
{
    return PropertyConditionList[i].GetPotentialFrequency(z, status);
};
double BoundaryConditions::GetPotentialInitialPhase(int i, double z, int& status) const
{
    return PropertyConditionList[i].GetPotentialInitialPhase(z, status);
};

std::vector<double> BoundaryConditions::GetPropertyConditionManualRestictions(int i)
{
    std::vector<double> result(3);
    result[0] = double(PropertyConditionList[i].AttachedElectrodeNumber);
    result[2] = double(PropertyConditionList[i].manualLineTag);
    result[1] = PropertyConditionList[i].manualLinePos;
    return result;
};
void BoundaryConditions::SetPropertyConditionManualRestictions(int i, std::vector<double> params)
{
    // PropertyConditionList[i].AttachedElectrodeNumber = int(params[0]);
    PropertyConditionList[i].manualLineTag = int(params[1]);
    PropertyConditionList[i].manualLinePos = params[0];
};
template bool
BoundaryConditions::CheckBoundaryConditionIntersection<double>(const DGeo::Edge<double>& edge,
                                                               int                       i);
template bool
BoundaryConditions::CheckBoundaryConditionIntersection<float>(const DGeo::Edge<float>& edge, int i);

template <class PointType>
bool BoundaryConditions::CheckBoundaryConditionIntersection(const DGeo::Edge<PointType>& edge,
                                                            int                          i)
{
    if (PropertyConditionList[i].manualLineTag == 0)
    {
        if (edge.point1.x > PropertyConditionList[i].manualLinePos &&
            edge.point2.x < PropertyConditionList[i].manualLinePos)
            return true;
        if (edge.point1.x < PropertyConditionList[i].manualLinePos &&
            edge.point2.x > PropertyConditionList[i].manualLinePos)
            return true;
    };
    if (PropertyConditionList[i].manualLineTag == 1)
    {
        if (edge.point1.y > PropertyConditionList[i].manualLinePos &&
            edge.point2.y < PropertyConditionList[i].manualLinePos)
            return true;
        if (edge.point1.y < PropertyConditionList[i].manualLinePos &&
            edge.point2.y > PropertyConditionList[i].manualLinePos)
            return true;
    };
    return false;
};

double PropertyCondition::GetPotentialOffset(double z, int& status) const
{
    int i;
    for (i = 0; i < zArray.size(); i++)
    {
        if (zArray[i][0] <= z && zArray[i][1] >= z)
            break;
    };
    if (i == zArray.size())
    {
        status = 0;
        return -1;
    };
    status = 1;
    return conditionProperties[i][0];
};

double PropertyCondition::GetPotentialAmplitude(double z, int& status) const
{
    int i;
    for (i = 0; i < zArray.size(); i++)
    {
        if (zArray[i][0] <= z && zArray[i][1] >= z)
            break;
    };
    if (i == zArray.size())
    {
        status = 0;
        return -1;
    };
    status = 1;
    return conditionProperties[i][1];
};
std::vector<std::string> PropertyCondition::GetPropertiesSimpleNames() const
{
    return propertiesNames;
};

double PropertyCondition::GetPotentialFrequency(double z, int& status) const
{
    int i;
    for (i = 0; i < zArray.size(); i++)
    {
        if (zArray[i][0] <= z && zArray[i][1] >= z)
            break;
    };
    if (i == zArray.size())
    {
        status = 0;
        return -1;
    };
    status = 1;
    return conditionProperties[i][2];
};
double PropertyCondition::GetPotentialInitialPhase(double z, int& status) const
{
    int i;
    for (i = 0; i < zArray.size(); i++)
    {
        if (zArray[i][0] <= z && zArray[i][1] >= z)
            break;
    };
    if (i == zArray.size())
    {
        status = 0;
        return -1;
    };
    status = 1;
    return conditionProperties[i][3];
};

void PropertyCondition::SetPropertiesFromFile(std::string file)
{
    FILE* fp = fopen(file.c_str(), "r");
    char  ss[250];

    int i = -1;
    int j = 1;
    zArray.clear();
    conditionProperties.clear();
    while (fgets(ss, 250, fp))
    {
        zArray.resize(j);
        zArray[j - 1].clear();
        conditionProperties.resize(j);
        double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
        int    tmpi;
        sscanf(ss, "%lf %lf %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6);
        zArray[j - 1].push_back(tmp1);
        zArray[j - 1].push_back(tmp2);
        conditionProperties[j - 1] = std::vector<double>{tmp3, tmp4, tmp5, tmp6};
        j++;
    }
    fclose(fp);
};

PropertyCondition::PropertyCondition(std::string typeIn, int boundaryTypeFlag)
{
    type = typeIn;
    conditionProperties.resize(1);
    conditionProperties[0]  = {0};
    Typeflag                = boundaryTypeFlag;
    manualLineTag           = 0;
    AttachedElectrodeNumber = 0;
    manualLinePos           = 0;
    if (typeIn == "potential")
        propertiesNames = fieldParameters;
    else
    {
        for (int i = 0; i < flowBoundaryTypeNames.size(); i++)
        {
            if (typeIn == flowBoundaryTypeNames[i])
            {
                propertiesNames = properiesNames[i];
                break;
            }
        };
    }
};
PropertyCondition::PropertyCondition()
{
    conditionProperties.resize(1);
    conditionProperties[0] = {0};
};

std::vector<std::vector<double>>& PropertyCondition::Get_zArray()
{
    return zArray;
};
std::vector<double> PropertyCondition::GetSimpleProperties() const
{
    return conditionProperties[0];
};

void PropertyCondition::SetSimpleProperties(std::vector<double> in)
{
    conditionProperties[0] = in;
};
double PropertyCondition::GetPotentialOffset() const
{
    return conditionProperties[0][0];
};

double PropertyCondition::GetPotentialAmplitude() const
{
    return conditionProperties[0][1];
};

double PropertyCondition::GetPotentialFrequency() const
{
    return conditionProperties[0][2];
};
double PropertyCondition::GetPotentialInitialPhase() const
{
    return conditionProperties[0][3];
};

double PropertyCondition::GetPotential(double t) const
{
    return conditionProperties[0][0] +
           conditionProperties[0][1] *
               std::cos(2 * commtools::PI() * conditionProperties[0][2] * t + conditionProperties[0][3]);
};
double PropertyCondition::GetPotential(double t, double z, int& status) const
{
    int i;
    for (i = 0; i < zArray.size(); i++)
    {
        if (zArray[i][0] <= z && zArray[i][1] >= z)
            break;
    };
    if (i == zArray.size())
    {
        status = 0;
        return -1;
    };
    status = 1;
    return conditionProperties[i][0] +
           conditionProperties[i][1] *
               std::cos(2 * commtools::PI() * conditionProperties[i][2] * t + conditionProperties[i][3]);
};

template void
PropertyCondition::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                         const unsigned int file_version) const;

template void
PropertyCondition::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                         const unsigned int file_version);

template <class Archive>
void PropertyCondition::save(Archive& ar, const unsigned int) const
{
    ar& conditionProperties;
    ar& boundariesList;
    ar& type;
    ar& Typeflag;
    ar& zArray;
    ar& manualLineTag;
    ar& AttachedElectrodeNumber;
    ar& manualLinePos;
    ar& propertiesNames;
}
template <class Archive>
void PropertyCondition::load(Archive& ar, const unsigned int)
{
    ar& conditionProperties;
    ar& boundariesList;
    ar& type;
    ar& Typeflag;
    ar& zArray;
    ar& manualLineTag;
    ar& AttachedElectrodeNumber;
    ar& manualLinePos;
    ar& propertiesNames;
}

int BoundaryConditions::PropertyConditionListSize() const
{
    return int(PropertyConditionList.size());
};
std::string BoundaryConditions::GetConditionPropertyType(int i) const
{
    return PropertyConditionList[i].type;
}
void BoundaryConditions::SetDefaultConditionsList(const std::vector<int>& in)
{
    DefaultConditionsList = in;
};
std::vector<int> BoundaryConditions::GetDefaultConditionsList() const
{
    return DefaultConditionsList;
};
std::vector<PropertyCondition> BoundaryConditions::GetPropertyConditionsList() const
{
    return PropertyConditionList;
};
void BoundaryConditions::AddDefaultConditionsList(int i)
{
    DefaultConditionsList.push_back(i);
};
void BoundaryConditions::SetPropertyConditionsBoundariesList(int i, const std::vector<int>& in)
{
    PropertyConditionList[i].boundariesList = in;
};
void BoundaryConditions::AddPropertyCondition(std::string type, int boundaryTypeFlag)
{
    PropertyConditionList.push_back(*(new PropertyCondition(type, boundaryTypeFlag)));
};
std::vector<int> BoundaryConditions::GetPropertyConditionsBoundariesList(int i) const
{
    return PropertyConditionList[i].boundariesList;
};
int BoundaryConditions::GetPropertyConditionTypeFlag(int i) const
{
    return PropertyConditionList[i].Typeflag;
};

int BoundaryConditions::GetNumberProperties() const
{
    return int(PropertyConditionList.size());
};

std::vector<double> BoundaryConditions::GetConditionPropertiesSimple(int i) const
{
    return PropertyConditionList[i].GetSimpleProperties();
};

std::vector<std::string> BoundaryConditions::GetConditionPropertiesSimpleNames(int i) const
{
    return PropertyConditionList[i].GetPropertiesSimpleNames();
};

void BoundaryConditions::SetConditionProperties(int i, std::vector<double> cond)
{
    PropertyConditionList[i].SetSimpleProperties(cond);
};
void BoundaryConditions::SetConditionPropertiesFromFile(int i, std::string cond)
{
    PropertyConditionList[i].SetPropertiesFromFile(cond);
};

template void
BoundaryConditions::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                          const unsigned int file_version) const;

template void
BoundaryConditions::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                          const unsigned int file_version);

template <class Archive>
void BoundaryConditions::save(Archive& ar, const unsigned int) const
{
    ar& DefaultConditionsList;
    ar& PropertyConditionList;
}
template <class Archive>
void BoundaryConditions::load(Archive& ar, const unsigned int)
{
    ar& DefaultConditionsList;
    ar& PropertyConditionList;
}