#ifndef PICA_BENCHMARK_UTILITY_PARTICLE_GENERATOR_H
#define PICA_BENCHMARK_UTILITY_PARTICLE_GENERATOR_H


#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleTraits.h"

#include <cmath>


namespace utility {

namespace detail {

template<class ParticleRef>
void generateParticle(ParticleRef particle, Random& random)
{
    particle.setMass(pica::constants::electronMass);
    particle.setCharge(pica::constants::electronCharge);
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
    double temperature = 1.0;
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


} // namespace utility::detail


template<class ParticleArray>
ParticleArray generateParticles(int numParticles)
{
    Random random;
    ParticleArray particles;
    for (int i = 0; i < numParticles; i++) {
        typename ParticleArray::Particle particle;
        detail::generateParticle(particle, random);
        particles.pushBack(particle);
    }
    return particles;
}

template<class Ensemble>
Ensemble generateParticles(pica::Int3 numCells, int numParticlesPerCell)
{
    typename Ensemble::PositionType minPosition(0.0, 0.0, 0.0);
    typename Ensemble::PositionType maxPosition(1.0, 1.0, 1.0);
    Ensemble particles(minPosition, maxPosition);
    Random random;
    long numParticles = numCells.volume() * numParticlesPerCell;
    for (int i = 0; i < numParticles; i++) {
        typename Ensemble::Particle particle;
        detail::generateParticle(particle, random);
        particles.add(particle);
    }
    return particles;
}


} // namespace utility


#endif
