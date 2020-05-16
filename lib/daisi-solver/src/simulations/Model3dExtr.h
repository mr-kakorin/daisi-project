#ifndef PIC3DEXTR_H
#define PIC3DEXTR_H
#include "DataTypes.h"

class Model3dExtrfloat : public ModelTemplate<device3dExtrfloat, float>
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
    Model3dExtrfloat(){
        //	simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    };
};
class Model3dExtrdouble : public ModelTemplate<device3dExtrdouble, double>
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
    Model3dExtrdouble(){
        //	simulationData.AddFlags(flagStringsSolver::simulationDataNamesBasePIC);
    };
};
#endif