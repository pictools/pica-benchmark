#include "pica/math/Constants.h"
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


template<class ParticleArray>
void process(ParticleArray& particles, int beginIdx, int endIdx)
{
    const Vector3<double> e(1.0, -1.0, 1.0);
    const Vector3<double> b(-1.0, 1.0, -1.0);
    const double dt = 0.1 / Constants<double>::c();
    BorisPusherBaseline<typename ParticleArray::Particle> pusher;
    #pragma simd
    #pragma forceinline
    for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
        pusher.push(particles[particleIdx], e, b, dt);
}

// Push all particles by one time step
template<class ParticleArray>
void runIteration(ParticleArray& particles)
{
    // Each thread processes some particles
    const int numParticles = particles.size();
    const int numThreads = pica::getNumThreads();
    const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
    #pragma omp parallel for
    for (int idx = 0; idx < numThreads; idx++) {
        const int beginIdx = idx * particlesPerThread;
        const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
        process(particles, beginIdx, endIdx);
    }
}

template<class ParticleArray>
void runBenchmark(ParticleArray& particles, int numIterations)
{
    for (int i = 0; i < numIterations; i++)
        runIteration(particles);
}

template<class ParticleArray>
ParticleArray generateParticles(int numParticles)
{
    ParticleArray result;
    for (int i = 0; i < numParticles; i++)
    {
        typename ParticleArray::Particle particle;
        result.pushBack(particle);
    }
    return result;
}

struct Parameters {
    int numParticles;
    int numIterations;
    int numThreads;
};

Parameters readParameters(int args, char* argv[])
{
    Parameters result;
    result.numParticles = 1000000;
    result.numIterations = 100;
    result.numThreads = 4;
    return result;
}

int main(int argc, char* argv[])
{
    Parameters parameters = readParameters(argc, argv);
    ParticleArrayAoS<ParticleAoS> particles =
        generateParticles<ParticleArrayAoS<ParticleAoS> >(parameters.numParticles);
    omp_set_num_threads(parameters.numThreads);
    runBenchmark(particles, parameters.numIterations);
    return 0;
}
