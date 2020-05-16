#ifndef DTLDevice_H
#define DTLDevice_H
#include "AccelDeviceBase.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

class DTLDevice : public AccelDeviceBase<AccelFlow>
{

  public:
    void AddFlow();
    int SetSomeParametersFromFile(int n, const std::string& filename, std::string& errorMessage, const std::string& folder);
    void TranslateParameters();
    void CreateSequences();
    DTLDevice();
};

#endif