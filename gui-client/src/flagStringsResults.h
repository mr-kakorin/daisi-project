#ifndef FLAGSTRINGSResults_H
#define FLAGSTRINGSResults_H
#include <string>
#include <vector>
namespace flagStrings
{
const static std::vector<std::string> PlotFlags2daxs = {"Er, V/m", "Ez, V/m", "Enorm, V/m",
                                                        "V, V",    "Bphi, T", "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags2d = {"Ex, V/m", "Ey, V/m", "Enorm, V/m", "V, V",
                                                     "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags2dpolar = {"Er, V/m", "Ephi, V/m", "Enorm, V/m", "V, V",
                                                          "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags3d = {
    "Ex, V/m", "Ey, V/m", "Ez, V/m", "Enorm, V/m", "Bx, V/m", "By, V/m", "Bz, V/m", "V, V", "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags3dPlane = {"YZ", "XZ", "XY"};

// const static std::vector<std::string> ResultsNames2daxs = { "RZ plane", "Energy (Z)", "Energy", "Energy phi", "XY
// plane", "RZ plane Rotating", "RZ plane Rotating 3 Trace", "R(t)", "Z(t)", "Phi(t)", "Charge(t)", "R(Z)" };
const static std::vector<std::string> ResultsNames2dpolar = {
    "XY plane", "X(t)", "Y(t)", "Currents", "Energy", "XY plane rotate", "Energy phi"};

const static std::vector<std::string> ResultsNames2d = {"Traces in XY plane", "X(t)",     "Y(t)", "PX(t)", "PY(t)",
                                                        "Charge(t)",          "Energy(t)"};
const static std::vector<std::string> ResultsNames2dLinac = {
    "Traces in XY plane", "X(t)", "Y(t)", "Z(t)",  "PX(t)", "PY(t)",
    "Energy(t)",          "X(Z)", "Y(Z)", "PX(Z)", "PY(Z)", "Energy(Z)"};

const static std::vector<std::string> ResultsNames2daxs = {
    "Traces in RZ plane", "Traces in ZR plane", "Traces in XY plane", "R(t)", "Z(t)", "PR(t)", "PZ(t)",
    "Charge(t)",          "Energy(t)"};

const static std::vector<std::string> ResultsNames2daxsLinac = {
    "Traces in RZ plane", "Traces in ZR plane", "R(t)", "Z(t)", "PR(t)", "PZ(t)", "Energy(t)", "R(Z)", "PR(Z)", "PZ(Z)",
    "Energy(Z)"};

const static std::vector<std::string> ResultsNames3d = {"X(t)", "Y(t)", "Z(t)", "X(z)", "Y(z)", "X(Y)", "Energy"};
};
#endif