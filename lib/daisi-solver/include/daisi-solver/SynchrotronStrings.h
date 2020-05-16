#ifndef SynchrotronStrings_H
#define SynchrotronStrings_H

#include <string>
#include <vector>

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

extern std::vector<std::string> SynchrotronTwissDynFlags;

extern std::vector<std::string> SynchrotronDynFlags;
extern std::vector<std::string> SynchrotronDynCMFlags;
extern std::vector<std::string> SynchrotronRecCMFlags;

extern std::vector<std::string> SynchrotronEstFlags;
extern std::vector<std::string> SynchrotronShuffleFlags;

extern std::string SynchrotronSolversNameCenterOfMass;
extern std::string SynchrotronSolversNameClosedOrbit;
extern std::string SynchrotronSolversNameTwiss;
extern std::string SynchrotronSolversBeam;
extern std::string SynchrotronSolversMADX;
extern std::string SynchrotronSolversCorrMatrix;
extern std::string SynchrotronSolversSVD;
extern std::string SynchrotronSolversMikado;
extern std::string SynchrotronSolversOrbitCalc;
extern std::string SynchrotronSolversAcceptance;
extern std::string SynchrotronSolversOptimization;
extern std::string SynchrotronSolversTolerancesNameReconstruction;
extern std::string SynchrotronSolversCorrEffName;
extern std::string SynchrotronSolversMagnetShuffleName;

extern std::vector<std::string>              SynchrotronOptixsflags;
extern std::vector<std::vector<std::string>> CorrectParametersNames;
extern std::vector<std::string>              SynchrotronCorrFlags;
extern std::vector<std::string>              SynchrotronOrbitFlags;
extern std::vector<std::string>              SynchrotronAccFlags;

#endif