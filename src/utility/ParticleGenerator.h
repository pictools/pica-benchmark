#ifndef PICA_BENCHMARK_UTILITY_PARTICLE_GENERATOR_H
#define PICA_BENCHMARK_UTILITY_PARTICLE_GENERATOR_H


#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"

#include <cmath>

namespace utility {

namespace detail
{
void initParticleTypes(int numParticleTypes)
{
    pica::ParticleType electron;
    electron.mass = pica::constants::electronMass;
    electron.charge = pica::constants::electronCharge;
    pica::ParticleType proton;
    proton.mass = pica::constants::protonMass;
    proton.charge = -pica::constants::electronCharge;
    // For now only electrons and protons are supported
    pica::ParticleTypes::typesVector.resize(numParticleTypes);
    for (int i = 0; i < numParticleTypes; i++)
        if (i % 2)
            pica::ParticleTypes::typesVector[i] = electron;
        else
            pica::ParticleTypes::typesVector[i] = proton;
    pica::ParticleTypes::types = &pica::ParticleTypes::typesVector[0];
    pica::ParticleTypes::numTypes = numParticleTypes;
}

template<class Particle>
void generateParticle(Particle& particle, Random& random, int numParticleTypes)
{
    particle.setType(rand() % numParticleTypes);
    particle.setFactor(1.0);
    typedef pica::Vector3<double> PositionType; /// todo - remove hardcoded type
    PositionType position;
    PositionType minPosition(0.0, 0.0, 0.0);
    PositionType maxPosition(1.0, 1.0, 1.0);
    for (int d = 0; d < pica::VectorDimensionHelper<PositionType>::dimension; d++)
        position[d] = minPosition[d] + (maxPosition[d] - minPosition[d]) * random.getUniform();
    particle.setPosition(position);
    // The standard deviation is sqrt(1/(2*alpha)), where alpha is
    // 3/2 * ((T/mc^2 + 1)^2 - 1)^(-1)
    double temperature = 1e-2 * pica::constants::electronMass * pica::constants::c * pica::constants::c;
    double alpha = temperature / particle.getMass() / pica::constants::c / pica::constants::c + 1;
    alpha = 1.5 / (alpha * alpha - 1);
    double sigma = sqrt(0.5 / alpha) * particle.getMass() * pica::constants::c;
    // Initial particle momentum is combination of given initial
    // momentum based on coords and random term in N(0, sigma)
    typedef pica::Vector3<double> MomentumType; /// todo - remove hardcoded type
    MomentumType momentum;
    for (int d = 0; d < pica::VectorDimensionHelper<MomentumType>::dimension; d++)
        momentum[d] = random.getNormal() * sigma;
    particle.setMomentum(momentum);
}

} // namespace detail

template<class ParticleArray>
ParticleArray generateParticles(int numParticles, int numParticleTypes)
{
    Random random;
    ParticleArray particles;
    detail::initParticleTypes(numParticleTypes);
    for (int i = 0; i < numParticles; i++) {
        typename ParticleArray::Particle particle;
        detail::generateParticle(particle, random, numParticleTypes);
        particles.pushBack(particle);
    }
    return particles;
}

template<class Ensemble>
Ensemble generateParticles(pica::Int3 numCells, int numParticlesPerCell, int numParticleTypes)
{
    typename Ensemble::PositionType minPosition(0.0, 0.0, 0.0);
    typename Ensemble::PositionType maxPosition(1.0, 1.0, 1.0);
    Ensemble particles(minPosition, maxPosition);
    Random random;
    long numParticles = numCells.volume() * numParticlesPerCell;
    detail::initParticleTypes(numParticleTypes);
    for (int i = 0; i < numParticles; i++) {
        typename Ensemble::Particle particle;
        detail::generateParticle(particle, random, numParticleTypes);
        particles.add(particle);
    }
    return particles;
}

template<class Ensemble>
Ensemble generateParticles(pica::Int3 numCells, int numParticlesPerCell, pica::Vector3<int> numSupercells, pica::Vector3<int> numSupercellsPerCell, int numParticleTypes)
{
    typename Ensemble::PositionType minPosition(0.0, 0.0, 0.0);
    typename Ensemble::PositionType maxPosition(1.0, 1.0, 1.0);
    Ensemble particles(minPosition, maxPosition, numSupercells, numSupercellsPerCell);
    Random random;
    long numParticles = numCells.volume() * numParticlesPerCell;
    detail::initParticleTypes(numParticleTypes);
    for (int i = 0; i < numParticles; i++) {
        typename Ensemble::Particle particle;
        detail::generateParticle(particle, random, numParticleTypes);
        particles.add(particle);
    }
    return particles;
}

} // namespace utility


#endif
