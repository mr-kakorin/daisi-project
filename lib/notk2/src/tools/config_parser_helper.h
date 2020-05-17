#ifndef NOTK_CONFIG_PARSER_HELPER
#define NOTK_CONFIG_PARSER_HELPER

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <memory>

#define CHECKVALUE(op, result, str)                                                                \
    if (op)                                                                                        \
    {                                                                                              \
        str    = str + #op + "; ";                                                                 \
        result = false;                                                                            \
    }

#endif
