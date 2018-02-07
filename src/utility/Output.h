#ifndef PICA_BENCHMARK_UTILITY_OUTPUT_H
#define PICA_BENCHMARK_UTILITY_OUTPUT_H


#include <string>


namespace utility {


struct PusherParameters;

void printHeader(const std::string& message, const PusherParameters& parameters);
void printResult(const PusherParameters& parameters, double runTimeSeconds);


} // namespace utility


#endif
