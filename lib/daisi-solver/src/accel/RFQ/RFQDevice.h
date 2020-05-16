#ifndef RFQDevice_H
#define RFQDevice_H

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "../base/AccelDeviceBase.h"
#include "../base/AccelFlow.h"

class RFQDevice : public AccelDeviceBase<AccelFlow>
{

  public:
    void AddFlow();
    int SetSomeParametersFromFile(int n, const std::string& filename, std::string& errorMessage,
                                  const std::string& folder);
    void TranslateParameters();
    void CreateSequences();
    RFQDevice();
    bool checkSequence();
};

#endif