#ifndef RFQModel_H
#define RFQModel_H

#include <AccelModelInterface.h>

#include "../base/AccelModelTemplate.h"

class RFQDevice;
class RFQSolver;
class RFQModel : public AccelModelTemplate<RFQDevice, RFQSolver>
{

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    RFQModel(const std::string& dataFolderIn)
        : AccelModelTemplate<RFQDevice, RFQSolver>(dataFolderIn){};
    void Solve(const std::string& solverName, double& progress, bool& flagAbort,
               std::string& errorMsg, std::vector<std::string>& status);
    void GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                    std::vector<std::string>&         names);
};

#endif