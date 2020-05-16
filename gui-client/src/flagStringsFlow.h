#ifndef FLAGSTRINGSFLOW_H
#define FLAGSTRINGSFLOW_H
#include <string>
#include <vector>
namespace flagStrings
{
const static std::vector<std::string> NamesOfNumbersParams2daxs = {"Emit each timestep", "At each position",
                                                                   "Along emitter"};
const static std::vector<std::string> NamesOfDistribParamsParams2daxs = {"Plasma temperature, eV",
                                                                         "Plasma density, cm^-3"};

const static std::vector<std::string> NamesOfNumbersParamsLinac       = {"Particles per bunch", "Number of bunches"};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac = {"Average bunch current, A",      // 0
                                                                         "Average energy, eV ",           // 1
                                                                         "Momentum spread, %",            // 2
                                                                         "Bunches frequency, Hz",         // 3
                                                                         "Min phase, rad",                // 4
                                                                         "Max phase, rad",                // 5
                                                                         "EmittanceX (norm), pi*mm*mrad", // 6
                                                                         "X max, m",                      // 7
                                                                         "dX/dZ max, rad",                // 8
                                                                         "EmittanceY (norm), pi*mm*mrad", // 9
                                                                         "Y max, m",                      // 10
                                                                         "dY/dZ max, rad",                // 11
                                                                         "Center mass pos., X, m",        // 12
                                                                         "Center mass pos., Y, m",        // 13
                                                                         "Center mass pos., dX/dZ, rad",  // 14
                                                                         "Center mass pos., dY/dZ, rad"}; // 15

const static std::vector<std::string> NamesOfNumbersParams2d       = {};
const static std::vector<std::string> NamesOfDistribParamsParams2d = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac2d       = {};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac2d = {};

const static std::vector<std::string> NamesOfNumbersParams2dpolar       = {};
const static std::vector<std::string> NamesOfDistribParamsParams2dpolar = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac2dpolar       = {};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac2dpolar = {};

const static std::vector<std::string> NamesOfNumbersParams3dfromExtrusion       = {};
const static std::vector<std::string> NamesOfDistribParamsParams3dfromExtrusion = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac3dfromExtrusion = {
    "Emission type", "Emit period", "Energy distribution", "Particles per step"};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac3dfromExtrusion = {
    "Average bunch current, A", "Average energy, eV ", "Beam radius, m"};

const static std::vector<std::string> radioDistributionStyleNames = {"Emission from fixed plasma surface",
                                                                     "Evaporation from irradiated material",
                                                                     "Field emission",
                                                                     "Thermionic emission",
                                                                     "Photoemission",
                                                                     "Accelerator-style injection"};
const static std::vector<std::string> radioEmissionTypeNames = {"Space-charge limited", "Uniform and user defined"};

const static std::vector<std::string> radioParticleTypeNames = {"Electon", "Ion"};
};
#endif