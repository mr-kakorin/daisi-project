#ifndef Dtools_H
#define Dtools_H
#include <vector>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

void strsplit(char* string, char* split, std::vector<double>& result, std::string& error);
void strsplit(char* string, char* split, std::vector<std::string>& result);
void strsplit(char* string, char* split, std::vector<float>& result, std::string& error);
void strsplit(std::string string, char* split, std::vector<std::string>& result);
class myunsorted_map
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& keys;
        ar& values;
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& keys;
        ar& values;
    }
    std::vector<std::string> keys;
    std::vector<double>      values;

  public:
	void clear();
    std::vector<std::string>& GetKeys();
    std::vector<double>&      GetValues();
    void SetValues(const std::vector<double>& in);
    void insert(const std::string& v1, double v2);
    int find(const std::string& key, double& val);
	void set(const std::string& key, const double val);

    double find(const std::string& key);
    //~myunsorted_map();
};


void GetTime(std::string& time, std::string& date);

template <typename T>
std::vector<T> json_vector(const boost::property_tree::ptree& pt, const std::string& key)
{
	std::vector<T> r;
	for (auto& item : pt.get_child(key))
	{
		r.push_back(item.second.get_value<T>());
	}
	return r;
}
bool read_property_tree(const boost::property_tree::ptree& config,
	const std::string&                 error_msg_prefix,
	const std::function<bool(const boost::property_tree::ptree&)>& operation);

std::shared_ptr<boost::property_tree::ptree> readJSONString(const std::string& config);

std::shared_ptr<boost::property_tree::ptree> readJSONFile(const std::string& path);

#endif