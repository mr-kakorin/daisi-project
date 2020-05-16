#ifndef BOOSTLOGCONF_H
#define BOOSTLOGCONF_H

#include <ostream>
#include <string>
#include <vector>

#include <boostlog0/boostlog.h>

std::ostream& operator<<(std::ostream& strm, sev_lvl level);

namespace initboostlog
{
void init(const sev_lvl sevLvl, const std::string&                generalFilePostfix,
          const std::vector<std::pair<std::string, std::string>>& channels_files = {},
          const bool exclude_from_general = true);
}

#endif // BOOSTLOGCONF_H
