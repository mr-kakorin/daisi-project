#ifndef BOUNDARYCONDITIONSGPU_H
#define BOUNDARYCONDITIONSGPU_H
#include "BoundaryConditions.h"
#include "GPUTypes.h"
#include <cuda_runtime.h>
class PropertyConditionGPU
{
  public:
    std::vector<double>   conditionProperties;
    GPUTypes::vector<int> boundariesList;
    GPUTypes::string      type;
    __device__ __host__ PropertyConditionGPU(){

    };
    __device__ __host__ PropertyConditionGPU(PropertyCondition& PropertyConditionCPU){
        //	conditionProperties[0] = PropertyConditionCPU.conditionProperties[0];
        //	boundariesList = PropertyConditionCPU.boundariesList;
        //	type = GPUTypes::string(PropertyConditionCPU.type);
    };
    /*	PropertyCondition(std::string typeIn)
            {
                    type = typeIn;
                    conditionProperties[0] = 0;
            };
            PropertyCondition()
            {
                    conditionProperties[0] = 0;
            };*/
};
class BoundaryConditionsGPU
{
    GPUTypes::vector<int>                  DefaultConditionsList;
    GPUTypes::vector<PropertyConditionGPU> PropertyConditionList;

  public:
    BoundaryConditionsGPU(){

    };
    __device__ __host__ BoundaryConditionsGPU(BoundaryConditions& BoundaryConditionsCPU)
    {
        DefaultConditionsList = BoundaryConditionsCPU.GetDefaultConditionsList();
        PropertyConditionList = BoundaryConditionsCPU.GetPropertyConditionsList();
    };
    /*int PropertyConditionListSize()
    {
            return int(PropertyConditionList.size());
    };
    void SetconditionProperties[0](int i, double cond)
    {
            PropertyConditionList[i].conditionProperties[0] = cond;
    }
    std::string GetconditionProperties[0]Type(int i)
    {
            return PropertyConditionList[i].type;
    }
    double GetconditionProperties[0](int i)
    {
            return PropertyConditionList[i].conditionProperties[0];
    }
    void SetDefaultConditionsList(const std::vector<int>& in)
    {
            DefaultConditionsList = in;
    };
    std::vector<int>& GetDefaultConditionsList()
    {
            return DefaultConditionsList;
    };
    void AddDefaultConditionsList(int i)
    {
            DefaultConditionsList.push_back(i);
    };
    void SetPropertyConditionsBoundariesList(int i, const std::vector<int> & in)
    {
            PropertyConditionList[i].boundariesList = in;
    };
    void AddPropertyCondition(std::string type)
    {
            PropertyConditionList.push_back(*(new PropertyCondition(type)));
    };
    std::vector<int> GetPropertyConditionsBoundariesList(int i)
    {
            return PropertyConditionList[i].boundariesList;
    };
    int GetNumberProperties()
    {
            return int(PropertyConditionList.size());
    };*/
};

#endif