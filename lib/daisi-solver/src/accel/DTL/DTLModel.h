#ifndef DTLModel_H
#define DTLModel_H

#include <AccelModelInterface.h>

#include "../base/AccelModelTemplate.h"

class DTLDevice;
class DTLSolver;
class DTLModel : public AccelModelTemplate<DTLDevice, DTLSolver>
{

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    DTLModel(const std::string& dataFolderIn)
        : AccelModelTemplate<DTLDevice, DTLSolver>(dataFolderIn){};
    void Solve(const std::string& solverName, double& progress, bool& flagAbort,
               std::string& errorMsg, std::vector<std::string>& status);
    void GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                    std::vector<std::string>&         names);
};

#endif