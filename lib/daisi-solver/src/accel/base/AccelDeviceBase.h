#ifndef AccelDeviceBase_H
#define AccelDeviceBase_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
class AccelFlow;
class myunsorted_map;
template <class FlowType>
class AccelDeviceBase
{

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const;
    template <class Archive> void load(Archive& ar, const unsigned int);

  protected:
    std::shared_ptr<myunsorted_map>              mainAccelParameters;
    std::shared_ptr<myunsorted_map>              mainAccelParametersFlags;
    std::vector<std::shared_ptr<myunsorted_map>> SectionParameters;
    std::vector<std::shared_ptr<FlowType>>      flows;
    std::vector<std::string>                     filesNames;
    std::vector<std::string>                     filesExtensions;
    std::vector<std::string>                     files;
    std::vector<std::string>                     namesOfSequences;
    std::vector<std::vector<double>>             Sequences;
    std::vector<std::string>                     calcParametersNames;
    std::vector<std::vector<double>>             calcParameters;
    std::shared_ptr<myunsorted_map>              calcMainParameters;

  public:
    std::vector<std::string> GetcalcMainParameters();
    std::vector<double>      GetMainAccelParameterCalculated();
    void checkInputParameters(std::string& errorMessage);
    std::vector<std::string>& GetSomeFileName();
    void GetBrouseFlags(std::vector<std::string>& brouse, std::vector<std::string>& names);
    std::vector<std::shared_ptr<FlowType>>& GetLinacFlows();
    void GetMainAccelParameters(std::vector<std::string>& keysO, std::vector<double>& p);
    void GetMainAccelParametersFlags(std::vector<std::string>& keysO, std::vector<double>& p);
    std::vector<double>& GetMainAccelParameters();
    void SetMainAccelParameters(const std::vector<double>& in);
    double GetParameter(const std::string& key);
    void GetSectionParameters(std::vector<std::vector<std::string>>& keysO, std::vector<std::vector<double>>& p);
    void SetAccelSectionParameters(int currentAccelSection, const std::vector<double>& params);
    std::vector<std::string>& GetnamesOfSequences();
    void GetSequencesOfParameters(int seq, std::vector<double>& sequence);
    void GetSequencesOfParametersCalc(int seq, std::vector<double>& sequence);
    std::vector<std::vector<double>>& GetSequencesOfParameters();
    std::vector<std::vector<double>>& GetSequencesOfParametersCalc();
    std::vector<std::string>&         GetcalcParametersNames();

    //	~AccelDeviceBase();
};
#endif
