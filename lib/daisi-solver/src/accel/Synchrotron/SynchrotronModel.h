#ifndef SynchrotronModel_H
#define SynchrotronModel_H

#include <AccelModelInterface.h>

#include "../base/AccelModelTemplate.h"


class SynchrotronDevice;
class SynchrotronSolver;
class SynchrotronModel : public AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>
{

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const;
    template <class Archive> void load(Archive& ar, const unsigned int);

  public:
    SynchrotronModel(const std::string& dataFolderIn)
        : AccelModelTemplate<SynchrotronDevice, SynchrotronSolver>(dataFolderIn){};
    void Solve(const std::string& solverName, double& progress, bool& flagAbort,
               std::string& errorMsg, std::vector<std::string>& status);
    void GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                    std::vector<std::string>&         names);
};

#endif