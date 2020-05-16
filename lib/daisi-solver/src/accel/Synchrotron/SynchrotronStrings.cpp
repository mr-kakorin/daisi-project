#include "SynchrotronStrings.h"

// std::vector <std::vector <  std::string >> SynchrotronInputFiles = { {},{},{},{ "" },{},{},{},{}
// }; std::vector <std::vector <  std::string >> SynchrotronInputFilesProps = { {},{},{} ,{ "*.madx"
// },{},{},{},{} }; std::vector <std::vector <  std::string >> SynchrotronInputFilesExt = { {},{},{}
// ,{ "*.madx" },{},{},{},{} }; std::vector <std::vector <  std::string >>
// SynchrotronInputFilesDescr = { {},{},{} ,{ "MAD-X config file"
// },{},{},{},{} };

// std::vector <std::vector <  std::string >> SynchrotronOutputFiles = { {},{ "twiss_daisi.txt" },{
// "trackone_daisi.txt"
// },{ "twiss_madx.txt","trackone_madx.txt"  },{"cororbitSVD.dat"},{ "cororbitMikado.dat" },{},{} };
// std::vector <std::vector <  std::string >> SynchrotronOutputFilesDescr = { {},{ "file with twiss
// parameters" }, { "file with dynamics" },{ "file with twiss parameters","file with dynamics" },{
// "file with corrections results" },{ "file with corrections results" },{},{} };

// std::vector <std::vector <  std::string >> SynchrotronSolverParameteresNames = { {},{},{ "Number
// of visualized traces", "Number of saved traces to file" },{},{}, { "Number of correctors
// (X)","Number of correctors (Y)" },{
//},{"Number of particles", "dY/dZ seach range", "dX/dZ seach range", "Relative channel radius" } };

std::vector<std::string> SynchrotronTwissDynFlags = {"BETA",

                                                     "ALFA", "MU"};

std::vector<std::string> SynchrotronDynFlags = {"X(L)",
                                                "dX/dZ(L)",
                                                "Y(L)",
                                                "dY/dZ(L)",
                                                "Distribution in plane XdX",
                                                "Distribution in plane YdY",
                                                "Distribution in plane XY",
                                                "Transmission",
                                                "BETA",
                                                "ALFA",
                                                "Xcm and Ycm",
                                                "dXcm/dZ and dYcm/dZ"};
std::vector<std::string> SynchrotronDynCMFlags = {"X and Y(L)", "dX/dZ and dY/dZ (L)"};
std::vector<std::string> SynchrotronRecCMFlags = {"X(L)",    "dX/dZ(L)",  "Y(L)",    "dY/dZ(L)",
                                                  "Fitness", "BETA X",    "Alpha X", "BETA Y",
                                                  "Alpha Y", "Parameters"};

std::vector<std::string> SynchrotronEstFlags = {
    "BETX max. no corr", "BETY max. no corr", "X max. no corr",  "Y max. no corr",
    "BETX max. corr",    "BETY max. corr",    "X max. corr",     "Y max. corr",
    "Xcorr/X",           "Ycorr/Y",           "X max. corr, lt", "Y max. corr, lt"};

std::vector<std::string> SynchrotronShuffleFlags = {"X max.", "Y max." };

std::string SynchrotronSolversNameCenterOfMass             = "Dynamics Sim, center of mass";
std::string SynchrotronSolversNameTwiss                    = "Dynamics Sim, Twiss";
std::string SynchrotronSolversBeam                         = "Dynamics Sim, beam";
std::string SynchrotronSolversMADX                         = "Dynamics Sim, MADX";
std::string SynchrotronSolversCorrMatrix                   = "Correction matrix calculation";
std::string SynchrotronSolversSVD                          = "Orbit correction (SVD)";
std::string SynchrotronSolversMikado                       = "Orbit correction (Mikado)";
std::string SynchrotronSolversOrbitCalc                    = "Orbit calculation";
std::string SynchrotronSolversAcceptance                   = "Acceptance calculation";
std::string SynchrotronSolversOptimization                 = "Orbit optimization";
std::string SynchrotronSolversTolerancesNameReconstruction = "Reconstruction of tolerances";
std::string SynchrotronSolversCorrEffName                  = "Estimation correction efficiency";
std::string SynchrotronSolversNameClosedOrbit = "Closed orbit calculation";
std::string SynchrotronSolversMagnetShuffleName = "Magnets shuffle";

std::vector<std::string> SynchrotronOptixsflags = {"DRIFT", "RBEND",     "QUADRUPOLE", "KICKER",
                                                   "SBEND", "SEXTUPOLE", "MARKER",     "MONITOR"};
std::vector<std::vector<std::string>> CorrectParametersNames = {
    {"L"}, {"L", "ANGLE", "E1", "E2"}, {"L", "k1"}, {}, {}, {"L"}};
std::vector<std::string> SynchrotronCorrFlags = {"Currents for correction",
                                                 "X coordinate deviation", "Y coordinate deviation",
                                                 "Singular values"};
std::vector<std::string> SynchrotronOrbitFlags = {"Orbit X and Y"};
std::vector<std::string> SynchrotronAccFlags   = {"Accepted distribution in plane XdX",
                                                "Accepted distribution in plane YdY",
                                                "Accepted distribution in plane XY"};
