#ifndef PICA_BENCHMARK_UTILITY_OUTPUT_H
#define PICA_BENCHMARK_UTILITY_OUTPUT_H


#include <string>


namespace utility {


struct FullParameters;
struct PusherParameters;

void printHeader(const std::string& message, const FullParameters& parameters);
void printHeader(const std::string& message, const PusherParameters& parameters);

void printResult(const FullParameters& parameters, double runTimeSeconds);
void printResult(const PusherParameters& parameters, double runTimeSeconds);


} // namespace utility


#endif
