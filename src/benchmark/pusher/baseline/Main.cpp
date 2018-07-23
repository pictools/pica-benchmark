#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleArray.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>

using pica::FP;
using pica::FP3;

struct Particle
{
    FP3 position;
    FP3 momentum;
    int typeIndex;

    FP getMass() { return pica::ParticleTypes::types[typeIndex].mass; }
    FP getCharge() { return pica::ParticleTypes::types[typeIndex].charge; }
};

template<class ParticleArray, class FieldValue>
void runBenchmark(ParticleArray& particles,
    const FieldValue& electricFieldValue,
    const FieldValue& magneticFieldValue,
    const utility::PusherParameters& parameters);

int main(int argc, char* argv[])
{
    utility::PusherParameters parameters = utility::readPusherParameters(argc, argv);
    utility::printHeader("pusher-baseline benchmark: using baseline 3D Boris particle pusher implementation and AoS particle representation",
        parameters);

    // Generate particles randomly,
    // particular coordinates and other data are not important for this benchmark
    std::vector<Particle> particles(parameters.numParticles);
    utility::Random random;
    utility::detail::initParticleTypes(parameters.numParticleTypes);
    for (int i = 0; i < parameters.numParticles; i++) {
        pica::Particle3d particle;
        utility::detail::generateParticle(particle, random, parameters.numParticleTypes);
        particles[i].position = particle.getPosition();
        particles[i].momentum = particle.getMomentum();
        particles[i].typeIndex = particle.getType();
    }

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

    for (int i = 0; i < parameters.numIterations; i++) {
        // Each thread processes some particles
        const int numParticles = particles.size();
        const int numThreads = pica::getNumThreads();
        const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;

        #pragma omp parallel for
        for (int idx = 0; idx < numThreads; idx++) {
            const int beginIdx = idx * particlesPerThread;
            const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);

            #pragma omp simd
            #pragma forceinline
            for (int i = beginIdx; i < endIdx; i++) {
                FP3 eMomentum = electricFieldValue * particles[i].getCharge() * dt /
                    ((FP)2.0 * particles[i].getMass() * pica::Constants<FP>::c());
                FP3 um = particles[i].momentum / (particles[i].getMass() * pica::Constants<FP>::c()) + eMomentum;
                FP3 t = magneticFieldValue * particles[i].getCharge() * dt /
                    ((FP)2.0 * particles[i].getMass() * pica::Constants<FP>::c() * sqrt(1.0 + um.norm2()));
                FP3 uprime = um + cross(um, t);
                FP3 s = t * (FP)2.0 / ((FP)1.0 + t.norm2());
                particles[i].momentum = (um + pica::cross(uprime, s) + eMomentum) * particles[i].getMass() * pica::Constants<FP>::c();
                FP3 v = particles[i].momentum / sqrt(pica::sqr(particles[i].getMass()) + (particles[i].momentum / pica::Constants<FP>::c()).norm2());
                particles[i].position += v * dt;
            }
        }
    }
}
