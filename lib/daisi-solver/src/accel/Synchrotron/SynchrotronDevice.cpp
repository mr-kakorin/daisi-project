#include <memory>

#include "FlagStringsSolver.h"
#include "SynchrotronDevice.h"
#include "SynchrotronTools.h"
#include "Tools.h"
#include "string.h"

#include "../base/AccelFlow.h"

std::vector<std::string> SynchrotronMainParameters = {"Number of cicles",
                                                      "Orbir radius, m",
                                                      "field/current",
                                                      "B*rho, Tl*m",
                                                      "Frequency of betatron oscillations",
                                                      "X aperture, m",
                                                      "Y aperture, m",
                                                      "Relative error in optics"};
std::vector<std::string> FilesExtensionsNucl = {"Optics description, (*.opt)",
                                                "Sequence description, (*.line)"};

std::vector<std::string> tmp;

template void
SynchrotronDevice::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);
template void
SynchrotronDevice::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                              const unsigned int file_version);

template void
SynchrotronDevice::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                         const unsigned int file_version);

template void
SynchrotronDevice::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                         const unsigned int file_version) const;

template <class Archive>
void SynchrotronDevice::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelDeviceBase>(*this);
    ar& optics;
    ar& opticElementsSequence;
}
template <class Archive>
void SynchrotronDevice::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelDeviceBase>(*this);
    ar& optics;
    ar& opticElementsSequence;
    opticElementsSequence->mtrans();
    TranslateParameters();
}

void SynchrotronDevice::SaveMADXConfigFile(std::string file)
{
    FILE*  fp     = fopen(file.c_str(), "w");
    double energy = flows[0]->GetTotalEnergy() * 1e-9;
    double dp     = flows[0]->getMomentumSpread();

    std::vector<double> cm;
    flows[0]->GetMassCenterVector(cm);

    fprintf(fp, "beam,\nparticle = ion,\nenergy = %lf;\n", energy);

    fprintf(fp, "call, file = \"opt.opt\";\n");
    fprintf(fp, "call, file = \"line.line\";\n");
    fprintf(fp, "option, echo, RBARC = FALSE;\n");
    // fprintf(fp, "option, echo;\n");

    fprintf(fp, "use, period = %s;\n", opticElementsSequence->GetName().c_str());

    fprintf(fp, "EOPTION, seed = 123456789, add = true;\nESAVE, FILE = sbend_quads_errors.txt;\n");

    fprintf(fp, "\ncall, file = \"align.txt\";\ncall, file = \"beta0.txt\";\ntwiss, BETA0 = "
                "initial, save, file = "
                "twiss_mad.txt;\nptc_create_universe;\nptc_create_layout, model = 2, method = 6, "
                "nst = 10, exact = "
                "false;\nptc_twiss, file = "
                "twiss_ptc.txt;\n\n");

    fprintf(fp, "ptc_create_universe;\nptc_create_layout, model = 2, method = 6, "
                "nst = 10, exact = false; \n");

    fprintf(fp, "ptc_start, x = %lf, px = %lf, y = %lf, py = %lf;\n", cm[0], cm[1], cm[2], cm[3]);

    fprintf(fp, "call, file = \"obs.txt\";\nptc_track, icase = 5, ");

    fprintf(fp, "DELTAP = %lf, ", dp);

    fprintf(fp, "FILE = cm_, element_by_element, "
                "dump = true, onetable = true, turns = 1, ffile = 1, \n		norm_no = "
                "1;\nptc_track_end;\nptc_end;\n\n");

    fprintf(
        fp,
        "ptc_create_universe;\nptc_create_layout, model = 2, method = 6, "
        "nst = 10, exact = false; \n	call, file = \"init.txt\";\nptc_track, icase = 5, ");

    fprintf(fp, "DELTAP = %lf, ", dp);

    fprintf(fp, "FILE = beam_, element_by_element, "
                "dump = true, onetable = true, turns = 1, ffile = 1, \n		norm_no = "
                "1;\nptc_track_end;\nptc_end;\n\n");

    fclose(fp);
};

void SynchrotronDevice::InitSequenceWithErrors(){
    /**opticElementsSequenceErrors = *opticElementsSequence;
    opticElementsSequenceErrors->InsertErrors(GetParameter("Relative error in optics"));
    opticElementsSequenceErrors->mtrans(GetParameter("Orbir radius, m"), GetParameter("field /
    current"), GetParameter("B*rho"));*/
};

void SynchrotronDevice::AddFlow()
{
    flows.push_back(std::shared_ptr<SynchrotronFlow>(new SynchrotronFlow()));
};

void SynchrotronDevice::GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                                   std::vector<std::string>&         names)
{
    props.clear();
    names.clear();
    props.resize(3);
    for (int i = 0; i < GetOpticElementsSequence()->length(); i++)
    {
        if (GetOpticElementsSequence()->GetType(i) != "DRIFT")
        {
            props[0].push_back(GetOpticElementsSequence()->getParameter(i, "at"));
            props[1].push_back(GetOpticElementsSequence()->getParameter(i, "L"));
            names.push_back(GetOpticElementsSequence()->GetType(i));
            // props[2].push_back(device->GetOpticElementsSequence()->getParameter(i, "L"));
        };
    };
};

std::shared_ptr<OpticElementsSequence>& SynchrotronDevice::GetOpticElementsSequence()
{
    return opticElementsSequence;
};

void SynchrotronDevice::CreateSequences(){

};

SynchrotronDevice::SynchrotronDevice()
{
    filesNames      = {"Elements description", "Elements sequence"};
    filesExtensions = FilesExtensionsNucl;
    files.resize(2);
    mainAccelParameters   = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    opticElementsSequence = std::shared_ptr<OpticElementsSequence>(new OpticElementsSequence());
    // opticElementsSequenceErrors = std::shared_ptr<OpticElementsSequence>(new
    // OpticElementsSequence());

    for (int i = 0; i < SynchrotronMainParameters.size(); i++)
        mainAccelParameters->insert(SynchrotronMainParameters[i], 0.0);

    mainAccelParametersFlags = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    // mainAccelParametersFlags->insert("Add errors to optics", 0.0);

    calcMainParameters = std::shared_ptr<myunsorted_map>(new myunsorted_map());
    std::vector<std::string> calcMainParametersNames = {
        "Accelerator length", "Number of correctors", "Number of pick-up electrodes"};

    for (int i = 0; i < calcMainParametersNames.size(); i++)
        calcMainParameters->insert(calcMainParametersNames[i], 0.0);
};

void SynchrotronDevice::TranslateParameters()
{
    double tmp;
    mainAccelParameters->find("Number of cicles", tmp);
    numberofcicles = tmp;

    mainAccelParameters->find("Orbir radius, m", p);

    mainAccelParameters->find("field/current", lambdak);

    mainAccelParameters->find("B*rho, Tl*m", Brho);

    mainAccelParameters->find("Frequency of betatron oscillations", Qx);
};

int SynchrotronDevice::SetSomeParametersFromFile(int n, const std::string& filename,
                                                 std::string&       errorMessage,
                                                 const std::string& folder)
{
    if (n == 0)
        return loadOptics(filename, errorMessage);
    if (n == 1)
    {
        int res = loadSeq(filename, errorMessage);
        saveCorrectors(folder);
        saveErrors(folder);

        return res;
    }
};
void SynchrotronDevice::saveErrors(const std::string& folder)
{
    auto errors = opticElementsSequence->GetErrorsStruct();
    save_errors(errors, folder + "tolerances.json");
}

void SynchrotronDevice::saveCorrectors(const std::string& folder)
{
    std::ofstream f_out(folder + "monitors.dat");
    for (size_t i = 0; i < opticElementsSequence->length(); i++)
    {
        if (opticElementsSequence->GetType(i) == "MONITOR")
        {
            f_out << opticElementsSequence->GetLabel(i) << "\t" << 1 << "\t" << 1 << "\n";
        }
    }
    f_out.close();
}
int SynchrotronDevice::loadOptics(const std::string& fileName, std::string& errorMessage)
{
    /*__try
    {
    }
    __except (EXCEPTION_EXECUTE_HANDLER)
    {
            errorMessage = "Something wrong in opt file.";
            return 0;
    }*/
    files[0] = fileName;
    optics.clear();
    FILE* fp = fopen(fileName.c_str(), "r");
    char  ss[250];

    if (!fp)
    {
        errorMessage = "Unable to open file.";
        return 0;
        // this = NULL;
    };
    errorMessage.clear();

    int str = 0;
    while (fgets(ss, 250, fp))
    {
        if (strlen(ss) && ss[0] != '/' && ss[1] != '/' && ss[0] != '\n')
        {
            strsplit(ss, ": ,;=, ", tmp);
            optics.push_back(std::shared_ptr<OpticElement>(new OpticElement(tmp, errorMessage)));
            if (errorMessage.size())
            {
                optics.clear();
                errorMessage = "Error in string " + std::to_string(str) + ". " + errorMessage +
                               ". \nOptics are not loaded.";
                fclose(fp);
                return 0;
            };
        }
        str++;
    };

    /*errorMessage = "Something wrong in opt file.";
    return 0;*/

    fclose(fp);

    for (int i = 0; i < optics.size(); i++)
    {
        int k = findOpticsElem(optics, optics[i]->GetType());
        if (k != -1)
            optics[i]->copy(optics[k]);
    };
    return 1;
};
int SynchrotronDevice::loadSeq(const std::string& fileName, std::string& errorMessage)
{
    files[1] = fileName;
    FILE* fp = fopen(fileName.c_str(), "r");
    char  ss[250];
    if (!fp)
    {
        errorMessage = "Unable to open file.";
        return 0;
        // this = NULL;
    };
    std::vector<std::string> tmp;

    std::vector<std::shared_ptr<OpticElementsSequence>> tmpSeq;

    while (fgets(ss, 250, fp))
    {
        if (strlen(ss) && ss[0] != '/' && ss[1] != '/' && ss[0] != '\n')
        {
            strsplit(ss, ": ,;=, ()\t", tmp);
            if (tmp.size() > 1 && tmp[1] == "LINE") // ��������� ������ LINE
            {
                tmpSeq.push_back(
                    std::shared_ptr<OpticElementsSequence>(new OpticElementsSequence(tmp[0])));
                for (int j = 2; j < tmp.size(); j++)
                {
                    if (!tmpSeq.back()->insertelem(optics, tmpSeq, tmp[j]))
                    {
                        errorMessage =
                            "There is no element " + tmp[j] +
                            " inside *.opt file. \nOptic elements sequence is not created.";
                        fclose(fp);
                        return 0;
                    }
                };
            }
            else // ����������� ������ LINE
            {
                for (int j = 0; j < tmp.size(); j++)
                {
                    if (!tmpSeq.back()->insertelem(optics, tmpSeq, tmp[j]))
                    {
                        errorMessage =
                            "There is no element " + tmp[j] +
                            " inside *.opt file. \nOptic elements sequence is not created.";
                        fclose(fp);
                        return 0;
                    }
                }
            }
        };
    };

    for (int i = 0; i < tmpSeq.size(); i++)
    {
        double at = 0;
        for (int j = 0; j < tmpSeq[i]->length(); j++)
        {
            tmpSeq[i]->setParameter(j, "at", at);
            at = at + tmpSeq[i]->getParameter(j, "L");
        }
        tmpSeq[i]->SetL(at);
    };
    fclose(fp);
    opticElementsSequence->clear();

    *opticElementsSequence = *(tmpSeq.back());
    opticElementsSequence->SetFlags();

    std::vector<double> p = {opticElementsSequence->GetL(),
                             double(opticElementsSequence->findType("KICKER").size()),
                             double(opticElementsSequence->findType("MONITOR").size())};

    calcMainParameters->SetValues(p);

    opticElementsSequence->GenerateErrors(0);
    opticElementsSequence->mtrans();

    return 1;
}
