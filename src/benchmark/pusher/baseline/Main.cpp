#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleBaseline.h"
#include "pica/particles/ParticleArray.h"
#include "pica/particlePush/BorisPusherBaseline.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>


template<class ParticleArray, class FieldValue>
void runBenchmark(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    const utility::PusherParameters& parameters);

template<class ParticleArray, class FieldValue>
void runIteration(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    double dt);

template<class ParticleArray, class FieldValue>
void process(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    int beginIdx, int endIdx, double dt);


int main(int argc, char* argv[])
{
    utility::PusherParameters parameters = utility::readPusherParameters(argc, argv);
    utility::printHeader("pusher-baseline benchmark: using baseline 3D Boris particle pusher implementation and AoS particle representation",
        parameters);

    // Generate particles randomly,
    // particular coordinates and other data are not important for this benchmark
    typedef pica::ParticleBaseline<pica::Three> Particle;
    typedef pica::ParticleArrayAoS<Particle> Particles;
    Particles particles = utility::generateParticles<Particles>(parameters.numParticles);

    // Generate random field values for each particle
    // to ensure there is no compile-time substitution of fields,
    // particular values of field are not important for this benchmark
    typedef pica::Vector3<double> FieldValue;
    FieldValue electricFieldValues = utility::generateField<double>();
    FieldValue magneticFieldValues = utility::generateField<double>();

    std::auto_ptr<utility::Stopwatch> timer(utility::createStopwatch());
    timer->start();
    runBenchmark(particles, electricFieldValues, magneticFieldValues, parameters);
    timer->stop();

    utility::printResult(parameters, timer->getElapsed());

    return 0;
}


// Run the whole benchmark
template<class ParticleArray, class FieldValue>
void runBenchmark(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    const utility::PusherParameters& parameters)
{
    omp_set_num_threads(parameters.numThreads);

    // value of time step does not matter for this benchmark, so make it
    // (somewhat) random to guarantee no compiler substitution is done for it
    utility::Random random;
    const double dt = random.getUniform() / pica::Constants<double>::c();

    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(particles, electricFieldValue, magneticFieldValue, dt);
}


// Push all particles by one time step
template<class ParticleArray, class FieldValue>
void runIteration(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    double dt)
{
    // Each thread processes some particles
    const int numParticles = particles.size();
    const int numThreads = pica::getNumThreads();
    const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
    #pragma omp parallel for
    for (int idx = 0; idx < numThreads; idx++) {
        const int beginIdx = idx * particlesPerThread;
        const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
        process(particles, electricFieldValue, magneticFieldValue, beginIdx, endIdx, dt);
    }
}


// Process particles with indexes in [beginIdx, endIdx) range
template<class ParticleArray, class FieldValue>
void process(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    int beginIdx, int endIdx, double dt)
{
    pica::BorisPusherBaseline<typename ParticleArray::Particle> pusher;
    #pragma simd
    #pragma forceinline
    for (int i = beginIdx; i < endIdx; i++)
        pusher.push(&particles[i], electricFieldValue, magneticFieldValue, dt);
}
