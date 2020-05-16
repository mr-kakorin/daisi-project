#ifndef PIC2DAXS_H
#define PIC2DAXS_H
#include "DataTypes.h"
#include "ModelTemplate.h"
class Model2daxsfloat : public ModelTemplate<device2daxsfloat, float>
{

    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }

    Model2daxsfloat()
    {
        //		simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    }
};
class Model2daxsdouble : public ModelTemplate<device2daxsdouble,
                                              double> // ModelTemplate <device2daxsdouble,
                                                      // Solver<double>, device2daxsdoubleGPU>
{

    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    Model2daxsdouble()
    {
        //		simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    }
};

#endif