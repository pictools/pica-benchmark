#include "utility/Parser.h"

#include "utility/Parameters.h"

#include <cmdline/cmdline.h>
#include "pica/threading/OpenMPHelper.h"

#include <string>


namespace utility {


PusherParameters readPusherParameters(int argc, char* argv[])
{
    cmdline::parser parser;
    parser.add<int>("nparticles", 0, "number of particles", false, 10000);
    parser.add<int>("niterations", 0, "number of iterations", false, 100);
    if (pica::useOpenMP())
        parser.add<int>("nthreads", 0, "number of OpenMP threads, default value is based on system settings",
            false, pica::getNumThreads());
    parser.parse(argc, argv);
    PusherParameters parameters;
    parameters.numParticles = parser.get<int>("nparticles");
    parameters.numIterations = parser.get<int>("niterations");
    if (pica::useOpenMP())
        parameters.numThreads = parser.get<int>("nthreads");
    else
        parameters.numThreads = 1;
    return parameters;
}


} // namespace utility
