#include "utility/Output.h"

#include "utility/Parameters.h"

#include "pica/threading/OpenMPHelper.h"

#include <iostream>


namespace utility {


void printHeader(const std::string& message, const PusherParameters& parameters)
{
    if (!message.empty())
        std::cout << message << "\n";
    std::cout << "Parameters:\n";
    std::cout << "    Number of particles: " << parameters.numParticles << "\n";
    std::cout << "    Number of time iterations: " << parameters.numIterations << "\n";
    if (pica::useOpenMP())
        std::cout << "    " << parameters.numThreads << " OpenMP threads\n";
    else
        std::cout << "    OpenMP is disabled\n";
    std::cout << std::endl;
}


void printResult(const PusherParameters& parameters, double runTimeSeconds)
{
    std::cout << "Total run time: " << runTimeSeconds << " sec";
    long numParticleUpdates = parameters.numParticles * parameters.numIterations;
    double nsPerParticleUpdate = runTimeSeconds * 1e9 / numParticleUpdates;
    std::cout << ", " << nsPerParticleUpdate << " ns per particle update";
    std::cout << std::endl;
}


} // namespace utility
