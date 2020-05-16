#ifndef PIC2D_H
#define PIC2D_H
#include "DataTypes.h"

class Model2dfloat : public ModelTemplate<device2dfloat, float>
{

    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    Model2dfloat(){
        //	simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    };
};
class Model2ddouble : public ModelTemplate<device2ddouble, double>
{

    friend class boost::serialization::access;

  public:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ModelTemplate>(*this);
    }
    Model2ddouble(){
        //		simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    };
};

#endif