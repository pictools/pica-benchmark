#ifndef PICA_BENCHMARK_UTILITY_PARSER_H
#define PICA_BENCHMARK_UTILITY_PARSER_H


#include "utility/Parameters.h"


namespace utility {


PusherParameters readPusherParameters(int argc, char* argv[]);
FullParameters readFullParameters(int argc, char* argv[]);


} // namespace utility


#endif
