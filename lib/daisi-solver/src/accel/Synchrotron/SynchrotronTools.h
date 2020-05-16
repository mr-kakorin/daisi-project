#ifndef SynchrotronTypes_H
#define SynchrotronTypes_H
//#define ARMA_MAT_PREALLOC 25

#include "armadillo"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <map>

class BaseOpticElement
{

};

class OpticElement
{
    friend class OpticElementsSequence;
    friend bool isEqual_el(const OpticElement* A, const OpticElement* B);
    friend int findOpticsElem(const std::vector<OpticElement*>& optics, const std::string& elem);
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& label;
        ar& type;
        ar& parameters;
        ar& at;
        ar& alignment_errors;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& label;
        ar& type;
        ar& parameters;
        ar& at;
        ar& alignment_errors;
    }

    std::string label;
    std::string type;
    double      at;
    std::map<std::string, double> parameters;
    std::vector<double> alignment_errors;

  public:
    void GenerateErrors(const std::map<std::string, std::array<double, 8>>& errors,
                        std::default_random_engine&       generator,
                        std::normal_distribution<double>& distribution);
    void GenerateErrors(double RelError, std::default_random_engine& generator,
                        std::normal_distribution<double>& distribution);

    std::pair<std::vector<double>, std::vector<double>> GetErrorsBorders(double RelError) const;
    std::pair<std::vector<double>, std::vector<double>>
    GetErrorsBorders(const std::map<std::string, std::array<double, 8>>& errors) const;

    void InsertErrors(const std::string& flagError, const std::vector<double>& RelError,
                      volatile int i1, volatile int& i2);
    std::vector<double> GetErrors() const;
    void InsertErrors(const std::vector<double>& errors);
    std::string GetLabel() const
    {
        return label;
    };
    std::string GetType() const
    {
        return type;
    };
	void SetLabel(const std::string& labe_l)
	{
		label = labe_l;
	}

    OpticElement(const std::vector<std::string>& input, std::string& error);

    OpticElement(const std::string type_)
        : type(type_){

          };

    void copy(const OpticElement& obj);
    void copy(const std::shared_ptr<OpticElement> obj);
    OpticElement* copy() const;
    OpticElement();
};
class OpticElementsSequence
{
    friend int findline(const std::vector<OpticElementsSequence>& tmpSeq, const std::string& str);
    friend class boost::serialization::access;
    //	friend void SynchrotronSolver::mtrans(const OpticElementsSequence& seq, double rorbit,
    // double lambdak, double
    // Brho, std::vector<arma::mat>& Tx, std::vector<arma::mat>& Ty);
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& elements;
        ar& name;
        ar& L;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& elements;
        ar& name;
        ar& L;
    }
    std::vector<OpticElement>        elements;
    double                           L;
    std::string                      name;
    std::default_random_engine       generator_loc;
    std::normal_distribution<double> distribution_loc;

  public:
	void insert_shuffle(const std::string& type, const std::vector<std::string>& shuffled, const std::map<std::string, double>& actual_vals);


    void SaveMADAlignmentCommands(const std::string& name);

    template <class T>
    void GenerateErrors(const T& errors)
    {
        for (int i = 0; i < length(); i++)
        {
            elements[i].GenerateErrors(errors, generator_loc, distribution_loc);
        }
    }

    std::vector<std::pair<std::string, std::array<double, 8>>> GetErrorsStruct() const;

    void InsertErrors(const std::vector<std::pair<std::string, std::array<double, 8>>>& errors);

    std::pair<size_t, size_t> GetErrIndexes(const std::vector<int>& monitiors) const;

    std::vector<double> GetErrors() const;

    template <class T>
    std::pair<std::vector<double>, std::vector<double>> GetErrorsBorders(const T& error) const
    {
        std::pair<std::vector<double>, std::vector<double>> result;
        for (int i = 0; i < length(); i++)
        {
            auto tmp = elements[i].GetErrorsBorders(error);
            result.first.insert(result.first.end(), tmp.first.begin(), tmp.first.end());
            result.second.insert(result.second.end(), tmp.second.begin(), tmp.second.end());
        }

        return result;
    }

    std::vector<arma::mat> Tx;
    std::vector<arma::mat> Ty;
    std::vector<arma::vec> Bx;
    std::vector<arma::vec> By;

    std::vector<arma::mat> sTx;
    std::vector<arma::mat> sTy;
    std::vector<arma::mat> sBx;
    std::vector<arma::mat> sBy;

    std::vector<arma::mat> sT;
    std::vector<arma::mat> T;
    std::vector<arma::mat> sB;
    std::vector<arma::mat> B;

    void clear()
    {
        elements.clear();
    };
    std::string GetName()
    {
        return name;
    };
    void SaveMADObsCommands(const std::string& name);
    void insertMonitors(int nmonitors);
    void SetFlags();
    void SetL(double Lin);
    double GetL() const;
    int setParameter(int i, std::string marker, double value);
    double getParameter(int i, std::string marker) const;
    size_t      length() const;
    std::string GetType(int i) const;
    std::string GetLabel(int i) const;
	void SetLabel(const int i, const std::string& label);

    OpticElementsSequence(const std::string& name);
    OpticElementsSequence();
    std::vector<int> findType(const std::string& name);
	std::vector<std::string> get_types_array(const std::string& type);
    void mtrans();
    int insertelem(std::vector<std::shared_ptr<OpticElement>>&          elementsBase,
                   std::vector<std::shared_ptr<OpticElementsSequence>>& tmpSeq,
                   const std::string&                                   name);
};

class LinacDynamicsAccel;
class SimulationDataAccel;
class SynchrotronDevice;
int findOpticsElem(std::vector<std::shared_ptr<OpticElement>>& optics, const std::string& elem);
void setOutputData(int i, int circle, const std::vector<int>& ind, const std::vector<int>& indP,
                   std::vector<void*>                                                       dataIn,
                   std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& data);
int findline(std::vector<std::shared_ptr<OpticElementsSequence>>& tmpSeq, const std::string& str);
arma::vec twiss(const arma::mat& M, const arma::mat& x);
void updateBeamPositions(const arma::mat& Mx, const arma::mat& My, LinacDynamicsAccel& dyn);
void savetwiss(const std::string& fileName, const std::vector<arma::mat>& xtwiss,
               const std::vector<arma::mat>& ytwiss, const std::string& particleType,
               double GeVmass, std::shared_ptr<OpticElementsSequence> seq,
               const std::vector<float>& mu_x, const std::vector<float>& mu_y);
void savetrack(const std::string& fileName, std::shared_ptr<OpticElementsSequence> seq,
               const std::vector<int>&                                           ind,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& madXData);
void loadtrack(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData);
void loadtwiss(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData);
void loadcm(const std::string& fileName, std::shared_ptr<SimulationDataAccel>& outputData);
void InitOrbit(std::shared_ptr<SynchrotronDevice>& device, std::vector<arma::mat>& x,
               std::vector<arma::mat>& y);
void loadtwiss(const std::string& fileName, float& betaM, float& alphaM);
void save_errors(const std::vector<std::pair<std::string, std::array<double, 8>>>& errors,
                 const std::string& fileName);
std::shared_ptr<std::vector<std::pair<std::string, std::array<double, 8>>>>
read_errors(const boost::property_tree::ptree& pt);
#endif
