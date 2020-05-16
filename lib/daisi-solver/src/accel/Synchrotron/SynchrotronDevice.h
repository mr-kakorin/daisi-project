#ifndef SynchrotronDevice_H
#define SynchrotronDevice_H

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "../base/AccelDeviceBase.h"
#include "../base/AccelFlow.h"

class OpticElement;
class OpticElementsSequence;
class SynchrotronDevice : public AccelDeviceBase<SynchrotronFlow>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template <class Archive> void save(Archive& ar, const unsigned int) const;
    template <class Archive> void load(Archive& ar, const unsigned int);

    int                                        numberofcicles;
    double                                     p, lambdak, Brho, Qx;
    std::vector<std::shared_ptr<OpticElement>> optics;
    std::shared_ptr<OpticElementsSequence>     opticElementsSequence;
    // std::shared_ptr<OpticElementsSequence> opticElementsSequenceErrors;
    int findelem(const std::string& elem) const;

  public:
    void                                    InitSequenceWithErrors();
    void                                    AddFlow();
    std::shared_ptr<OpticElementsSequence>& GetOpticElementsSequence();
    SynchrotronDevice();
    void saveCorrectors(const std::string& folder);
    void saveErrors(const std::string& folder);

    int loadOptics(const std::string& fileName, std::string& errorMessage);
    int loadSeq(const std::string& fileName, std::string& errorMessage);
    int SetSomeParametersFromFile(int n, const std::string& filename, std::string& errorMessage,
                                  const std::string& folder);
    void TranslateParameters();
    void CreateSequences();
    void GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                    std::vector<std::string>&         names);
    void SaveMADXConfigFile(std::string file);
};

#endif