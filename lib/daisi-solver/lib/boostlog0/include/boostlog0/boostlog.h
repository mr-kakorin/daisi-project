#ifndef BOOSTLOG0_H
#define BOOSTLOG0_H

#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/trivial.hpp>

#define BL_VAL(str) #str
#define BL_TOSTRING(str) BL_VAL(str)

enum class sev_lvl
{
    trace,
    debug,
    info,
    warning,
    error,
    critical
};

typedef boost::log::sources::severity_channel_logger_mt<sev_lvl, std::string> sev_ch_log;

BOOST_LOG_INLINE_GLOBAL_LOGGER_INIT(logger, sev_ch_log)
{
    return sev_ch_log();
}

#ifdef LOG_CHANNEL
#define LOG(lvl) BOOST_LOG_CHANNEL_SEV(logger::get(), BL_TOSTRING(LOG_CHANNEL), lvl)
#endif

#ifndef LOG_CHANNEL
#define LOG(lvl) BOOST_LOG_CHANNEL_SEV(logger::get(), "Undef.Chan.", lvl)
#endif

#define LOG_CH(lvl, chan) BOOST_LOG_CHANNEL_SEV(logger::get(), chan, lvl)

#define BL_TRACE() LOG(sev_lvl::trace)
#define BL_DEBUG() LOG(sev_lvl::debug)
#define BL_INFO() LOG(sev_lvl::info)
#define BL_WARNING() LOG(sev_lvl::warning)
#define BL_ERROR() LOG(sev_lvl::error)
#define BL_CRITICAL() LOG(sev_lvl::critical)

#define BL_TRACE_CH(chan) LOG_CH(sev_lvl::trace, chan)
#define BL_DEBUG_CH(chan) LOG_CH(sev_lvl::debug, chan)
#define BL_INFO_CH(chan) LOG_CH(sev_lvl::info, chan)
#define BL_WARNING_CH(chan) LOG_CH(sev_lvl::warning, chan)
#define BL_ERROR_CH(chan) LOG_CH(sev_lvl::error, chan)
#define BL_CRITICAL_CH(chan) LOG_CH(sev_lvl::critical, chan)

#endif
