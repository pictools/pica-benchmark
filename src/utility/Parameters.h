#ifndef PICA_BENCHMARK_UTILITY_PARAMETERS_H
#define PICA_BENCHMARK_UTILITY_PARAMETERS_H


namespace utility {


// Parameters for particle pusher benchmarks
struct PusherParameters {
    int numParticles;
    int numIterations;
    int numThreads;
};


} // namespace utility


#endif
