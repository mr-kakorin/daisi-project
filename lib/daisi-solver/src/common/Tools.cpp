#include "Tools.h"
#include <time.h>

std::string tmpStr;


void strsplit(std::string string, char* split, std::vector<std::string>& result)
{
	result.clear();
	char* pch = strtok(&string[0], split);

	while (pch != NULL)
	{
		result.push_back(pch);
		pch = strtok(NULL, split);
	}
}

void GetTime(std::string& timeIn, std::string& dateIn)
{
    time_t rawtime;
    tm*    timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    timeIn.clear();
    timeIn = "\"";
    dateIn.clear();
    dateIn = "\"";

    if (timeinfo->tm_mday < 10)
        dateIn = dateIn + "0";

    dateIn = dateIn + std::to_string(timeinfo->tm_mday);

    dateIn = dateIn + "/";

    if (timeinfo->tm_mon + 1 < 9 < 10)
        dateIn = dateIn + "0";

    dateIn = dateIn + std::to_string(timeinfo->tm_mon + 1);

    dateIn = dateIn + "/";

    dateIn = dateIn + std::to_string(timeinfo->tm_year - 100);

    dateIn = dateIn + "\"";

    if (timeinfo->tm_hour < 10)
        timeIn = timeIn + "0";

    timeIn = timeIn + std::to_string(timeinfo->tm_hour);

    timeIn = timeIn + ".";

    if (timeinfo->tm_min < 10)
        timeIn = timeIn + "0";

    timeIn = timeIn + std::to_string(timeinfo->tm_min);

    timeIn = timeIn + ".";

    if (timeinfo->tm_sec < 10)
        timeIn = timeIn + "0";

    timeIn = timeIn + std::to_string(timeinfo->tm_sec);

    timeIn = timeIn + "\"";

    /*if (timeinfo->tm_mon + 1<9)
            fprintf(fp, "@ DATE             %4s \"%d/0%d/%d\" \n", "%08s", timeinfo->tm_mday, timeinfo->tm_mon + 1,
    timeinfo->tm_year - 100); else fprintf(fp, "@ DATE             %4s \"%d/%d/%d\" \n", "%08s", timeinfo->tm_mday,
    timeinfo->tm_mon + 1, timeinfo->tm_year - 100);

    if (timeinfo->tm_min>9)
            fprintf(fp, "@ TIME             %4s \"%d\.%d\.%d\" \n", "%08s", timeinfo->tm_hour, timeinfo->tm_min,
    timeinfo->tm_sec); else fprintf(fp, "@ TIME             %4s \"%d\.0%d\.%d\" \n", "%08s", timeinfo->tm_hour,
    timeinfo->tm_min, timeinfo->tm_sec);*/
}

void strsplit(char* string, char* split, std::vector<std::string>& result)
{
    result.clear();
    char* pch = strtok(string, split);

    while (pch != NULL)
    {
        result.push_back(pch);
        pch = strtok(NULL, split);
    }
}

void strsplit(char* string, char* split, std::vector<double>& result, std::string& error)
{

    result.clear();
    char* pch = strtok(string, split);
    while (pch != NULL)
    {
        tmpStr = pch;
        result.push_back(std::stod(tmpStr));
        pch = strtok(NULL, split);
    }
}

void strsplit(char* string, char* split, std::vector<float>& result, std::string& error)
{

    result.clear();
    char* pch = strtok(string, split);
    while (pch != NULL)
    {
        tmpStr = pch;
        result.push_back(std::stof(tmpStr));
        pch = strtok(NULL, split);
    }
}

std::vector<std::string>& myunsorted_map::GetKeys()
{
    return keys;
}

std::vector<double>& myunsorted_map::GetValues()
{
    return values;
}

void myunsorted_map::SetValues(const std::vector<double>& in)
{
    values = in;
}

void myunsorted_map::insert(const std::string& v1, double v2)
{
    keys.push_back(v1);
    values.push_back(v2);
}

void myunsorted_map::clear()
{
	keys.clear();
	values.clear();
}

int myunsorted_map::find(const std::string& key, double& val)
{
    for (int i = 0; i < keys.size(); i++)
    {
        if (keys[i] == key)
        {
            val = values[i];
            return i;
        }
    }
    return -1;
}

void myunsorted_map::set(const std::string& key, const double val)
{
	for (int i = 0; i < keys.size(); i++)
	{
		if (keys[i] == key)
		{
			values[i] = val;
		}
	}
}


double myunsorted_map::find(const std::string& key)
{
    for (int i = 0; i < keys.size(); i++)
    {
        if (keys[i] == key)
            return values[i];
    }
    return -1;
}

bool read_property_tree(const boost::property_tree::ptree& config,
	const std::string&                 error_msg_prefix,
	const std::function<bool(const boost::property_tree::ptree&)>& operation)
{
	try
	{
		return operation(config);
	}

	catch (const boost::property_tree::json_parser_error& e)
	{
		BOOST_LOG_TRIVIAL(error) << error_msg_prefix << e.message();
		return false;
	}
	catch (const std::exception& e)
	{
		auto msg = e.what();
		BOOST_LOG_TRIVIAL(error) << error_msg_prefix << msg;
		return false;
	}
	catch (...)
	{
		BOOST_LOG_TRIVIAL(error) << error_msg_prefix << ": unknown error";
		return false;
	}
}

std::shared_ptr<boost::property_tree::ptree> readJSONString(const std::string& config)
{
	std::stringstream ss;
	ss << config;

	auto pt = std::make_shared<boost::property_tree::ptree>();
	try
	{
		boost::property_tree::read_json(ss, *pt);
		return pt;
	}
	catch (const boost::property_tree::json_parser_error& e)
	{
		BOOST_LOG_TRIVIAL(error) << "readJSONString::error: " << e.message();
		return nullptr;
	}
	catch (...)
	{
		BOOST_LOG_TRIVIAL(error) << "readJSONString::unknown error";
		return nullptr;
	}
}

std::shared_ptr<boost::property_tree::ptree> readJSONFile(const std::string& path)
{

	auto pt = std::make_shared<boost::property_tree::ptree>();
	try
	{
		boost::property_tree::read_json(path, *pt);
		return pt;
	}
	catch (const boost::property_tree::json_parser_error& e)
	{
		BOOST_LOG_TRIVIAL(error) << "readJSONFile::error: " << e.message();
		return nullptr;
	}
	catch (...)
	{
		BOOST_LOG_TRIVIAL(error) << "readJSONFile::unknown error";
		return nullptr;
	}
}