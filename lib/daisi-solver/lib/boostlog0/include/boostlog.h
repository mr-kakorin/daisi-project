#ifndef BOOSTLOG0_H
#define BOOSTLOG0_H

#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/trivial.hpp>

namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;

#define VAL(str) #str
#define TOSTRING(str) VAL(str)

enum class sev_lvl
{
    trace,
    debug,
    info,
    warning,
    error,
    critical
};

typedef src::severity_channel_logger_mt<sev_lvl, std::string> sev_ch_log;

BOOST_LOG_INLINE_GLOBAL_LOGGER_INIT(logger, sev_ch_log)
{
    return sev_ch_log();
}

#ifdef LOG_CHANNEL
#define LOG(lvl) BOOST_LOG_CHANNEL_SEV(logger::get(), TOSTRING(LOG_CHANNEL), lvl)
#endif

#ifndef LOG_CHANNEL
#define LOG(lvl) BOOST_LOG_CHANNEL_SEV(logger::get(), "UNDEF CHANNEL", lvl)
#endif

#define LOG_CH(lvl, chan) BOOST_LOG_CHANNEL_SEV(logger::get(), chan, lvl)

#endif