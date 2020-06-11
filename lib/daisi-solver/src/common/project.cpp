#include <iostream>

#include "project.h"
#ifdef SIM
#include "../simulations/Model2d.h"
#include "../simulations/Model2daxs.h"
#include "../simulations/Model2dpolar.h"
#include "../simulations/Model3dExtr.h"
#endif

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>

namespace Dproject {

    bool isCompatible(std::string currentVersion, std::string version) {
        return true;
    }

    void Dproject::project::SaveProject(std::string version) {
        std::string projectFilePath = (projectFolder + "/" + projectName);
        FILE *fpOut = fopen(projectFilePath.c_str(), "w");
        fprintf(fpOut, "%d\n%d\n%d\n%s\n%s\n", problemType, precisionType, solverType,
                currentModelFileName.c_str(), version.c_str());
        fclose(fpOut);
        SaveCurrentModel(currentModelFileName, version);
    }

    project::project(std::string inputFile, std::string version, std::string &errorMsg) {
        LoadFromFile(inputFile, version, errorMsg);
    }

    project::project(int problemTypeIn, int precisionTypeIn, int solverTypeIn,
                     std::string projectNameIn, std::string mainFolderIn, std::string version) {
        problemType = problemTypeIn;
        precisionType = precisionTypeIn;
        solverType = solverTypeIn;
        projectName = projectNameIn + ".dproj";
        version = version;
        projectFolder = mainFolderIn + "/" + projectNameIn;
        currentModelFileName = projectNameIn + "_0.mdl";

        meshParamFileName = projectFolder + "/mesh.dat";

        boost::filesystem::create_directory(projectFolder);
        boost::filesystem::create_directory(projectFolder + "/modelFiles");
        boost::filesystem::create_directory(projectFolder + "/dataFiles");

        FILE *fp = fopen(meshParamFileName.c_str(), "w");
        fclose(fp);

        currentModel = NULL;
        accelModel = NULL;
        ModelCreate();
        SaveProject(version);
    }

    void project::SaveData(std::string FileName) {
        FileName = (projectFolder + "/dataFiles" + "/" + FileName);
        //was nullptr
        accelModel->SaveData(FileName);
    }

    void project::LoadData(std::string FileName) {
        accelModel->LoadData(FileName);
    }

    void project::LoadProject(std::string FileName, std::string version, std::string &errorMsg) {
        this->version.assign(version);
        //	if (accelModel)
        //		delete accelModel;

        LoadFromFile(FileName, version, errorMsg);
    }

    void tryOpenProjectFile(std::string &inputFile, std::string &errorMsg) {
        FILE *fp = fopen(inputFile.c_str(), "r");
        char tmps[256];
        int tmp;
        if (!fp) {
            errorMsg = "Invalid file path.";
            return;
            // this = NULL;
        }
        fscanf(fp, "%d", &tmp);

        fscanf(fp, "%d", &tmp);

        fscanf(fp, "%d", &tmp);

        fscanf(fp, "%s", tmps, 256);

        fclose(fp);
    }

    void project::LoadFromFile(std::string inputFile, std::string version, std::string &errorMsg) {
        tryOpenProjectFile(inputFile, errorMsg);

        if (errorMsg.size())
            return;

        if (accelModel || currentModel)
            ModelDelete();

        FILE *fp = fopen(inputFile.c_str(), "r");

        if (!fp) {
            errorMsg = "Invalid file path.";
            return;
            // this = NULL;
        };

        char tmps[256];
        int tmp;

        fscanf(fp, "%d", &tmp);
        problemType = tmp;

        fscanf(fp, "%d", &tmp);
        precisionType = tmp;

        fscanf(fp, "%d", &tmp);
        solverType = tmp;

        fscanf(fp, "%s", tmps, 256);

        fclose(fp);

        currentModelFileName = tmps;

        size_t found = inputFile.find_last_of("//");

        projectFolder = inputFile.substr(0, found);

        meshParamFileName = projectFolder + "/mesh.dat";

        std::string modelFilePath = projectFolder + std::string("/modelFiles/") + currentModelFileName;
        projectName = inputFile.substr(found + 1);

        // currentModel = NULL;
        // accelModel = NULL;

        ModelCreate();

        ModelLoad(modelFilePath, version, errorMsg);
    }

    project::project() {
        currentModel = NULL;
        accelModel = NULL;
    }

    void project::SaveCurrentModel(std::string currentModelFileNameIn, std::string version) {
        std::string currentModelFilePath =
                (projectFolder + "/modelFiles" + "/" + currentModelFileNameIn);
        std::ofstream ofs(currentModelFilePath.c_str(), std::ios::out | std::ios::binary);
        boost::archive::binary_oarchive oa(ofs);

        oa << version;
        oa << problemType;
        oa << precisionType;

        switch (problemType) {

#ifdef SIM
            // case 8:
            // {
            //     DTLModel* ser1 = dynamic_cast<DTLModel*>(accelModel);
            //     oa << *ser1;
            // }
            // break;
            case 1:
                switch (precisionType) {
                    case 1: {
                        Model2dfloat *ser1 = dynamic_cast<Model2dfloat *>(currentModel);
                        oa << *ser1;
                    }
                        break;
                    case 0: {
                        Model2ddouble *ser2 = dynamic_cast<Model2ddouble *>(currentModel);
                        oa << *ser2;
                    }
                        break;
                }
                break;
            case 2:
                switch (precisionType) {
                    case 1: {
                        Model2daxsfloat *ser1 = dynamic_cast<Model2daxsfloat *>(currentModel);
                        oa << *ser1;
                    }
                        break;
                    case 0: {
                        Model2daxsdouble *ser2 = dynamic_cast<Model2daxsdouble *>(currentModel);
                        oa << *ser2;
                    }
                        break;
                }
                break;
            case 3:
                switch (solverType) {
                    case 0:
                        switch (precisionType) {
                            case 1: {
                                Model2dpolarfloat *ser1 = dynamic_cast<Model2dpolarfloat *>(currentModel);
                                oa << *ser1;
                            }
                                break;
                            case 0: {
                                Model2dpolardouble *ser2 = dynamic_cast<Model2dpolardouble *>(currentModel);
                                oa << *ser2;
                            }
                                break;
                        }
                        break;
                    case 1:
                        switch (precisionType) {
                            case 1: {
                                Model2dpolarfloat *ser1 = dynamic_cast<Model2dpolarfloat *>(currentModel);
                                oa << *ser1;
                            }
                                break;
                            case 0: {
                                Model2dpolardouble *ser2 = dynamic_cast<Model2dpolardouble *>(currentModel);
                                oa << *ser2;
                            }
                                break;
                        }
                        break;
                }
                break;

            case 4:
                switch (precisionType) {
                    case 1: {
                        Model3dExtrfloat *ser1 = dynamic_cast<Model3dExtrfloat *>(currentModel);
                        oa << *ser1;
                    }
                        break;
                    case 0: {
                        Model3dExtrdouble *ser2 = dynamic_cast<Model3dExtrdouble *>(currentModel);
                        oa << *ser2;
                    }
                        break;
                }
                break;
#endif
        }
    }

    void project::ModelLoad(std::string modelFilePath, std::string version, std::string &errorMsg) {

        size_t found = modelFilePath.find_last_of("/");

        currentModelFileName = modelFilePath.substr(found + 1);
        currentModelFilePath = modelFilePath;

        std::ifstream ifs(modelFilePath.c_str(), std::ios::in | std::ios::binary);
        try {
            boost::archive::binary_iarchive ia(ifs);

            std::string versionModel;
            int problemTypeModel;
            int precisionTypeModel;

            ia >> versionModel;
            ia >> problemTypeModel;
            ia >> precisionTypeModel;

            if (problemTypeModel != problemType || precisionTypeModel != precisionType) {
                errorMsg = "Incompatible project and model";
                return;
            }
            if (!isCompatible(version, versionModel)) {
                errorMsg = "Incompatible version of model file";
                return;
            }

            switch (problemType) {

#ifdef SIM

                // case 8:
                // {
                //     DTLModel* ser1 = dynamic_cast<DTLModel*>(accelModel);
                //     ia >> *ser1;
                // }
                // break;
                case 1:
                    switch (precisionType) {
                        case 1: {
                            Model2dfloat *ser1 = dynamic_cast<Model2dfloat *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                        case 0: {
                            Model2ddouble *ser1 = dynamic_cast<Model2ddouble *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                    }
                    break;
                case 2:
                    switch (precisionType) {
                        case 1: {
                            Model2daxsfloat *ser1 = dynamic_cast<Model2daxsfloat *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                        case 0: {
                            Model2daxsdouble *ser1 = dynamic_cast<Model2daxsdouble *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                    }
                    break;

                case 3:
                    switch (precisionType) {
                        case 1: {
                            Model2dpolarfloat *ser1 = dynamic_cast<Model2dpolarfloat *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                        case 0: {
                            Model2dpolardouble *ser1 = dynamic_cast<Model2dpolardouble *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                    }
                    break;
                case 4:
                    switch (precisionType) {
                        case 1: {
                            Model3dExtrfloat *ser1 = dynamic_cast<Model3dExtrfloat *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                        case 0: {
                            Model3dExtrdouble *ser1 = dynamic_cast<Model3dExtrdouble *>(currentModel);
                            ia >> *ser1;
                        }
                            break;
                    }
                    break;
#endif
            }
        }
        catch (std::exception &er) {
            std::string err = er.what();
            std::cout << er.what();
            int tt = 0;
        }
    }

    void project::ModelCreate() {

        switch (problemType) {
#ifdef SIM

            // case 8:
            //     accelModel = new DTLModel(projectFolder + "/dataFiles/");
            //     break;

            case 1:
                switch (precisionType) {
                    case 1: {
                        currentModel = new Model2dfloat();
                    }
                        break;
                    case 0: {
                        currentModel = new Model2ddouble();
                    }
                        break;
                }
                break;
            case 2:
                switch (precisionType) {
                    case 1: {
                        currentModel = new Model2daxsfloat();
                    }
                        break;
                    case 0: {
                        currentModel = new Model2daxsdouble();
                    }
                        break;
                }
                break;
            case 3:
                switch (precisionType) {
                    case 1: {
                        currentModel = new Model2dpolarfloat();
                    }
                        break;
                    case 0: {
                        currentModel = new Model2dpolardouble();
                    }
                        break;
                }
                break;
            case 4:
                switch (precisionType) {
                    case 1: {
                        currentModel = new Model3dExtrfloat();
                    }
                        break;
                    case 0: {
                        currentModel = new Model3dExtrdouble();
                    }
                        break;
                }
                break;
#endif
        };
    }

    void project::ModelDelete() {
        switch (problemType) {

#ifdef SIM
            // case 8:
            // {
            //     DTLModel* model = dynamic_cast<DTLModel*>(accelModel);
            //     delete model;
            // }
            // break;
            case 1:
                switch (precisionType) {
                    case 1: {
                        Model2dfloat *model = dynamic_cast<Model2dfloat *>(currentModel);
                        delete model;
                    }
                        break;
                    case 0: {
                        Model2ddouble *model = dynamic_cast<Model2ddouble *>(currentModel);
                        delete model;
                    }
                        break;
                }
                break;
            case 2:
                switch (precisionType) {
                    case 1: {
                        Model2daxsfloat *model = dynamic_cast<Model2daxsfloat *>(currentModel);
                        delete model;
                    }
                        break;
                    case 0: {
                        Model2daxsdouble *model = dynamic_cast<Model2daxsdouble *>(currentModel);
                        delete model;
                    }
                        break;
                }
                break;

            case 3:
                switch (precisionType) {
                    case 1: {
                        Model2dpolarfloat *model = dynamic_cast<Model2dpolarfloat *>(currentModel);
                        delete model;
                    }
                        break;
                    case 0: {
                        Model2dpolardouble *model = dynamic_cast<Model2dpolardouble *>(currentModel);
                        delete model;
                    }
                        break;
                }
                break;
            case 4:
                switch (precisionType) {
                    case 1: {
                        Model3dExtrfloat *model = dynamic_cast<Model3dExtrfloat *>(currentModel);
                        delete model;
                    }
                        break;
                    case 0: {
                        Model3dExtrdouble *model = dynamic_cast<Model3dExtrdouble *>(currentModel);
                        delete model;
                    }
                        break;
                }
                break;
#endif
        }
    }
}