#include <fstream>
#include <iostream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>

#include "boostlog.h"

// Define the attribute keywords
BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", sev_lvl)
BOOST_LOG_ATTRIBUTE_KEYWORD(channel, "Channel", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(object_attr, "Object", std::string)

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;

std::ostream &operator<<(std::ostream &strm, sev_lvl level)
{
    static const char *strings[] = {"trace", "debug", "info", "warning", "error", "critical"};

    if (static_cast<std::size_t>(level) < sizeof(strings) / sizeof(*strings))
    {
        strm << strings[static_cast<int>(level)];
    }
    else
    {
        strm << static_cast<int>(level);
    }

    return strm;
}
namespace initboostlog
{
    int  rotSize = 20 * 1024 * 1024;
    void addDataLogFile(const std::string &channelIn, sev_lvl sevLvl, const std::string &attr)
    {
        auto fmtDebug = expr::stream << expr::format_date_time<boost::posix_time::ptime>(
                                            "TimeStamp", "%Y-%m-%d %H:%M:%S")
                                     << ":\t" << expr::smessage;

        logging::add_file_log(keywords::file_name = "log/" + channelIn + "/" + attr + ".log",
                              keywords::rotation_size = rotSize, keywords::auto_flush = true,
                              keywords::format = fmtDebug);
    }

    void init(const sev_lvl sevLvl, const std::string &generalFilePostfix)
    {
        // формат общего логгера
        auto fmtGeneral =
            expr::stream
            << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
            << ": <" << expr::if_(expr::has_attr(severity))[expr::stream << severity]
            << expr::if_(!expr::has_attr(severity))[expr::stream << logging::trivial::severity]
            << "> " << expr::if_(expr::has_attr(channel))[expr::stream << "[" << channel << "] "]
            << expr::if_(!expr::has_attr(channel))[expr::stream << "["
                                                                << "Undef channel"
                                                                << "] "]
            << expr::smessage;

        // настройка тривиального лога

        auto lvlTriv = logging::trivial::trace;
        switch (sevLvl)
        {
            case sev_lvl::trace:
                lvlTriv = logging::trivial::trace;
                break;
            case sev_lvl::debug:
                lvlTriv = logging::trivial::debug;
                break;
            case sev_lvl::info:
                lvlTriv = logging::trivial::info;
                break;
            case sev_lvl::warning:
                lvlTriv = logging::trivial::warning;
                break;
            case sev_lvl::error:
                lvlTriv = logging::trivial::error;
                break;
            case sev_lvl::critical:
                lvlTriv = logging::trivial::fatal;
                break;
        }

        auto filter = ((!expr::has_attr(severity) && logging::trivial::severity >= lvlTriv) ||
                       (expr::has_attr(severity) && severity >= sevLvl));

        // файл с общим логом
        logging::add_file_log(keywords::file_name = "log/general_" + generalFilePostfix + ".log",
                              keywords::rotation_size = rotSize, keywords::filter = filter,
                              keywords::format = fmtGeneral, keywords::auto_flush = true);

        // консоль с общим логом
        logging::add_console_log(std::cout, keywords::format = fmtGeneral,
                                 keywords::filter = filter);

        logging::add_common_attributes();
    }
}