#ifndef BOOSTLOGCONF_H
#define BOOSTLOGCONF_H

#include <ostream>
#include <string>

#include <boostlog0/boostlog.h>

namespace channels
{
    const std::string System = "System";
}

std::ostream &operator<<(std::ostream &strm, sev_lvl level);

namespace initboostlog
{
    void init(const sev_lvl sevLvl, const std::string &generalFilePostfix);
}

#endif  // BOOSTLOGCONF_H