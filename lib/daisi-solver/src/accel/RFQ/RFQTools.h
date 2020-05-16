#ifndef RFQTOOLS_H
#define RFQTOOLS_H

#include "Dmath.h"

#include <common_tools/constants.h>

#include <boost/math/special_functions/bessel.hpp>
class AccelFlow;
int GenerateRFQForFlowPr(int& succes, const std::vector<double>& RFQParameters, std::shared_ptr<AccelFlow>& flow,
                         const std::vector<double>& Modulations, const std::vector<double>& SyncPhases,
                         const std::vector<double>& MatcherRadii, std::vector<double>& L,
                         std::vector<double>& MinumalRadii, std::vector<double>& MinumalRadiiRegular,
                         const std::vector<double>& OutRadii, std::vector<double>& AccEff,
                         std::vector<double>& AvEnergies);

#endif
