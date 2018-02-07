#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/math/Vectors.h"
#include "pica/particles/ParticleArray.h"
#include "pica/particlePush/BorisPusherBaseline.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>

using namespace pica;

// A simple type for particle
class ParticleAoS {
public:

    // Types for conforming ParticleInterface
    typedef double Real;
    typedef Vector3<Real> PositionType;
    typedef Vector3<Real> MomentumType;
    typedef Real GammaType;
    typedef Real MassType;
    typedef Real ChargeType;
    typedef float FactorType;

    ParticleAoS() :
        factor(1),
        mass(0.0),
        charge(0.0),
        invGamma(1.0) {}

    ParticleAoS(const PositionType& position, const MomentumType& momentum,
        MassType mass, ChargeType charge, FactorType factor = 1) :
        position(position), mass(mass), charge(charge), factor(factor)
    {
        setMomentum(momentum);
    }

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const
    {
        return p * Constants<MassType>::c() * mass;
    }

    void setMomentum(const MomentumType& newMomentum)
    {
        p = newMomentum / (Constants<GammaType>::c() * mass);
        invGamma = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
    }

    MomentumType getVelocity() const
    {
        return p * (Constants<GammaType>::c() * invGamma);
    }

    void setVelocity(const MomentumType& newVelocity)
    {
        p = newVelocity / sqrt(constants::c * constants::c - newVelocity.norm2());
        invGamma = (FP)1 / sqrt((FP)1 + p.norm2());
    }

    GammaType getGamma() const { return static_cast<GammaType>(1.0) / invGamma; }

    MassType getMass() const { return mass; }
    void setMass(MassType newMass) { mass = newMass; }

    ChargeType getCharge() const { return charge; }
    void setCharge(ChargeType newCharge) { charge = newCharge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

private:

    PositionType position;
    MomentumType p;
    MassType mass;
    ChargeType charge;
    FactorType factor;
    GammaType invGamma;
};



template<class ParticleArray, class FieldValue>
void runBenchmark(ParticleArray& particles,
    const std::vector<FieldValue>& electricFieldValues,
    const std::vector<FieldValue>& magneticFieldValues,
    const utility::PusherParameters& parameters);

template<class ParticleArray, class FieldValue>
void runIteration(ParticleArray& particles,
    const FieldValue* electricFieldValues,
    const FieldValue* magneticFieldValues,
    double dt);

template<class ParticleArray, class FieldValue>
void process(ParticleArray& particles,
    const FieldValue* electricFieldValues,
    const FieldValue* magneticFieldValues,
    int beginIdx, int endIdx, double dt);


int main(int argc, char* argv[])
{
    utility::PusherParameters parameters = utility::readPusherParameters(argc, argv);
    utility::printHeader("pusher-baseline benchmark: using baseline 3D Boris particle pusher implementation and SoA particle representation",
        parameters);

    typedef ParticleAoS Particle;
    typedef pica::ParticleArrayAoS<Particle> Particles;
    Particles particles = utility::generateParticles<Particles>(parameters.numParticles);

    typedef pica::Vector3<double> FieldValue;
    std::vector<FieldValue> electricFieldValues = utility::generateField<double>(parameters.numParticles);
    std::vector<FieldValue> magneticFieldValues = utility::generateField<double>(parameters.numParticles);

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
    const std::vector<FieldValue>& electricFieldValues,
    const std::vector<FieldValue>& magneticFieldValues,
    const utility::PusherParameters& parameters)
{
    omp_set_num_threads(parameters.numThreads);

    // value of time step does not matter for this benchmark, so make it
    // (somewhat) random to guarantee no compiler substitution is done for it
    utility::Random random;
    const double dt = random.getUniform() / pica::Constants<double>::c();

    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(particles, &electricFieldValues.front(),
            &magneticFieldValues.front(), dt);
}


// Push all particles by one time step
template<class ParticleArray, class FieldValue>
void runIteration(ParticleArray& particles,
    const FieldValue* electricFieldValues,
    const FieldValue* magneticFieldValues,
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
        process(particles, electricFieldValues, magneticFieldValues, beginIdx, endIdx, dt);
    }
}


// Process particles with indexes in [beginIdx, endIdx) range
template<class ParticleArray, class FieldValue>
void process(ParticleArray& particles,
    const FieldValue* electricFieldValues,
    const FieldValue* magneticFieldValues,
    int beginIdx, int endIdx, double dt)
{
    BorisPusherBaseline<typename ParticleArray::Particle> pusher;
    #pragma simd
    #pragma forceinline
    for (int i = beginIdx; i < endIdx; i++)
        pusher.push(particles[i], electricFieldValues[i], magneticFieldValues[i], dt);
}
