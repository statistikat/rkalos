#ifndef KALOS_VERSION_HPP
#define KALOS_VERSION_HPP
#if __has_include("kalos/revision.hpp")
#include "kalos/revision.hpp"
#endif
#ifndef KALOS_CORE_REVISION
#define KALOS_CORE_REVISION "development"
#endif
namespace kalos {
inline constexpr const char* version = "0.1.0";
inline constexpr const char* core_revision = KALOS_CORE_REVISION;
}
#endif
