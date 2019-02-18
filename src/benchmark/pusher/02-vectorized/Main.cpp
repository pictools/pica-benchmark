#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleArray.h"
#include "pica/particlePush/BorisPusher.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>

using pica::FP;
using pica::FP3;

struct Particle
{
    FP3 position;
    FP3 p;
    int typeIndex;
};

template<class ParticleArray, class FieldValue>
void runBenchmark(ParticleArray& particles,
    const FieldValue& electricFieldValues,
    const FieldValue& magneticFieldValues,
    const utility::PusherParameters& parameters);

int main(int argc, char* argv[])
{
    utility::PusherParameters parameters = utility::readPusherParameters(argc, argv);
    utility::printHeader("pusher-precomputing benchmark: using optimized 3D Boris particle pusher implementation with precomputing of inverse gamma and AoS particle representation",
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
        particles[i].p = particle.getP();
        particles[i].typeIndex = particle.getType();
    }

    // Generate random field values for each particle
    // to ensure there is no compile-time substitution of fields,
    // particular values of field are not important for this benchmark
    typedef FP3 FieldValue;
    FieldValue electricFieldValue = utility::generateField<double>();
    FieldValue magneticFieldValue = utility::generateField<double>();

    std::auto_ptr<utility::Stopwatch> timer(utility::createStopwatch());
    timer->start();
    runBenchmark(particles, electricFieldValue, magneticFieldValue, parameters);
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

    std::vector<double> coeff(parameters.numParticleTypes);
    for (int i = 0; i < parameters.numParticleTypes; i++)
    {
        coeff[i] = pica::ParticleTypes::typesVector[i].charge * dt /
            (2.0 * pica::ParticleTypes::typesVector[i].mass * pica::Constants<FP>::c());
    }

    for (int i = 0; i < parameters.numIterations; i++) {
        // Each thread processes some particles
        const int numParticles = particles.size();
        const int numThreads = pica::getNumThreads();
        const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;

        #pragma omp parallel for
        for (int idx = 0; idx < numThreads; idx++) {
            const int beginIdx = idx * particlesPerThread;
            const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);

            #pragma simd
            #pragma forceinline
            for (int i = beginIdx; i < endIdx; i++) {
                // The code below uses precomputed coefficient:
                // eCoeff = q * dt / (2 * m * c)
                FP eCoeff = coeff[particles[i].typeIndex];
                FP3 eMomentum = electricFieldValue * eCoeff;
                FP3 um = particles[i].p + eMomentum;
                FP3 t = magneticFieldValue * (eCoeff / sqrt((double)1.0 + um.norm2()));
                FP3 uprime = um + pica::cross(um, t);
                FP3 s = t * ((double)2.0 / ((double)1.0 + t.norm2()));
                particles[i].p = um + pica::cross(uprime, s) + eMomentum;
                FP3 v = particles[i].p / sqrt(1.0 + particles[i].p.norm2());
                particles[i].position += v * dt;
            }
        }
    }
}
