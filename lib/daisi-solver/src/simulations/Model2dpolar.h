#ifndef PIC2DPOLAR_H
#define PIC2DPOLAR_H
#include "ModelTemplate.h"

class Model2dpolarfloat : public ModelTemplate<device2dpolarfloat, float>
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

    Model2dpolarfloat()
    {
        //	simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    }
};
class Model2dpolardouble : public ModelTemplate<device2dpolardouble, double>
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

    Model2dpolardouble()
    {
        // simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    }
};

#endif