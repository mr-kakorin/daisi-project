#include <fstream>
#include <iostream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/attributes/value_extraction.hpp>
#include <boost/log/attributes/value_extraction_fwd.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/char_decorator.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/phoenix/bind/bind_function.hpp>

#include <boostlogconf/boostlog.h>

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
    std::vector<std::string> strings{"trace", "debug", "info", "warning", "error", "critical"};

    if (static_cast<size_t>(level) < strings.size())
    {
        strm << "<" << strings[static_cast<int>(level)] << "> ";
        for (size_t i = 0; i < 7 - strings[static_cast<int>(level)].size(); i++)
        {
            strm << " ";
        }
    }
    else
    {
        strm << static_cast<int>(level);
    }

    return strm;
}

/* expr::stream &operator<<(expr::stream &strm, logging::trivial::severity_level level)
{
} */
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
    std::string format_channel(logging::value_ref<std::string> const &filename)
    {
        std::stringstream out;
        out << "[" << filename.get() << "] ";
        for (size_t i = 0; i < 11 - filename.get().size(); i++)
        {
            out << " ";
        }
        return out.str();
    }

    bool check_channels(logging::value_ref<std::string> const &channel_name,
                        const std::vector<std::string> &       channels)
    {
        for (const auto &ch : channels)
        {
            if (channel_name.get() == ch)
            {
                return true;
            }
        }
        return false;
    }
    bool check_channel(logging::value_ref<std::string> const &channel_name,
                       const std::string &                    channel)
    {
        return channel_name.get() == channel;
    }
    void init(const sev_lvl sevLvl, const std::string &               generalFilePostfix,
              const std::vector<std::pair<std::string, std::string>> &channels_files)
    {
        // формат общего логгера
        auto fmtGeneral =
            expr::stream
            << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
            << ": " << expr::if_(expr::has_attr(severity))[expr::stream << severity]
            << expr::if_(!expr::has_attr(
                   severity))[expr::stream << "<" << logging::trivial::severity << "> "]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::trace)[expr::stream << std::setw(3)]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::debug)[expr::stream << std::setw(3)]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::info)[expr::stream << std::setw(4)]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::warning)[expr::stream << std::setw(1)]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::error)[expr::stream << std::setw(3)]
            << expr::if_(!expr::has_attr(severity) &&
                         logging::trivial::severity ==
                             logging::trivial::fatal)[expr::stream << std::setw(3)]
            << expr::if_(expr::has_attr(channel))[expr::stream << boost::phoenix::bind(
                                                      &format_channel,
                                                      expr::attr<std::string>("Channel"))]
            << expr::if_(!expr::has_attr(channel))[expr::stream << "["
                                                                << "Undef.Chan."
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

        std::vector<std::string> exclude_ch;
        for (const auto &ch : channels_files)
        {
            exclude_ch.push_back(ch.first);
        }

        auto filter = ((!expr::has_attr(severity) && logging::trivial::severity >= lvlTriv) ||
                       (expr::has_attr(severity) && severity >= sevLvl)) &&
                      (expr::has_attr(channel) &&
                       !boost::phoenix::bind(&check_channels, expr::attr<std::string>("Channel"),
                                             exclude_ch));

        // файл с общим логом
        logging::add_file_log(keywords::file_name = "log/general_" + generalFilePostfix + ".log",
                              keywords::rotation_size = rotSize, keywords::filter = filter,
                              keywords::format = fmtGeneral, keywords::auto_flush = true);

        // консоль с общим логом
        logging::add_console_log(std::cout, keywords::format = fmtGeneral,
                                 keywords::filter = filter);

        for (const auto &ch : channels_files)
        {
            auto filter_loc =
                ((!expr::has_attr(severity) && logging::trivial::severity >= lvlTriv) ||
                 (expr::has_attr(severity) && severity >= sevLvl)) &&
                (expr::has_attr(channel) &&
                 boost::phoenix::bind(&check_channel, expr::attr<std::string>("Channel"),
                                      ch.first));

            logging::add_file_log(keywords::file_name = "log/" + ch.second + ".log",
                                  keywords::rotation_size = rotSize, keywords::filter = filter_loc,
                                  keywords::format = fmtGeneral, keywords::auto_flush = true);
        }

        logging::add_common_attributes();
    }
}