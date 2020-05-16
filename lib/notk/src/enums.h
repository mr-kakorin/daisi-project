#ifndef NOTK_ENUMS_H
#define NOTK_ENUMS_H

#include <serreflection/defines.hpp>

namespace notk
{
enum class AccuracyType
{
    RELATIVE,
    ABSOLUTE,
    ZERO
};
enum class AccuracySourceType
{
    ARGUMENT,
    FITNESS
};
enum class InitialBordersType
{
    MANUAL,
    VECTOR
};
}
SERIALIZIBLE_ENUM(notk::AccuracyType, (RELATIVE)(ABSOLUTE)(ZERO))
SERIALIZIBLE_ENUM(notk::AccuracySourceType, (ARGUMENT)(FITNESS))
SERIALIZIBLE_ENUM(notk::InitialBordersType, (MANUAL)(VECTOR))

#endif
