#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/fieldSolver/YeeSolver.h"
#include "pica/grid/YeeGrid.h"
#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleArray.h"
#include "pica/particlePush/BorisPusherBaseline.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>


template<class Ensemble, class Grid>
void runBenchmark(Ensemble& particles, Grid& fields, 
    const utility::FullParameters& parameters);

template<class Ensemble, class Grid>
void runIteration(Ensemble& particles, Grid& fields, double dt);


int main(int argc, char* argv[])
{
    utility::FullParameters parameters = utility::readFullParameters(argc, argv);
    utility::printHeader("full-baseline benchmark: using unordered 3D particle ensemble, CIC form factor and SoA particle representation",
        parameters);

    // Generate particles randomly,
    // particular coordinates are not important for this benchmark
    typedef pica::Particle<pica::Three> Particle;
    typedef pica::ParticleArraySoA<Particle> ParticleArray;
    typedef typename pica::Ensemble<ParticleArray,
        pica::EnsembleRepresentation_Unordered>::Type Particles;
    Particles particles = utility::generateParticles<Particles>(
        parameters.numCells, parameters.particlesPerCell);

    // Generate fields
    typedef pica::YeeGrid<pica::Three> Grid;
    Grid fields = utility::generateField<Grid>(particles.getMinPosition(),
        particles.getMaxPosition(), parameters.numCells);

    std::auto_ptr<utility::Stopwatch> timer(utility::createStopwatch());
    timer->start();
    runBenchmark(particles, fields, parameters);
    timer->stop();

    utility::printResult(parameters, timer->getElapsed());

    return 0;
}


// Run the whole benchmark
template<class Ensemble, class Grid>
void runBenchmark(Ensemble& particles, Grid& fields,
    const utility::FullParameters& parameters)
{
    omp_set_num_threads(parameters.numThreads);

    // time step
    const double dt = 0.1 / pica::Constants<double>::c();

    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(particles, fields, dt);
}


// Simulate one time step
template<class Ensemble, class Grid>
void runIteration(Ensemble& particles, Grid& fields, double dt)
{
    updateParticles(particles, fields, dt);
    updateField(fields, dt);
}

template<class Ensemble, class Grid>
void updateParticles(Ensemble& particles, const Grid& fields, double dt)
{
    // Each thread processes some particles
    const int numParticles = particles.size();
    const int numThreads = pica::getNumThreads();
    const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
    #pragma omp parallel for
    for (int idx = 0; idx < numThreads; idx++) {
        const int beginIdx = idx * particlesPerThread;
        const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
        process(particles, fields, beginIdx, endIdx, dt);
    }
}

template<class Ensemble, class Grid>
void process(Ensemble& particles, const Grid& fields,
    int beginIdx, int endIdx, double dt)
{
    push(particles, fields, beginIdx, endIdx, dt);
    applyBoundaryConditions(particles, beginIdx, endIdx);
}

template<class Ensemble, class Grid>
void push(Ensemble& particles, const Grid& fields,
    int beginIdx, int endIdx, double dt)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle> particleArray = particles.getParticles();
    pica::BorisPusherBaseline<Particle> pusher;
    pica::FieldInterpolatorCIC<Grid> fieldInterpolator(fields);
    #pragma simd
    #pragma forceinline
    for (int i = beginIdx; i < endIdx; i++) {
        pica::Vector3<double> e, b;
        fieldInterpolator.get(particleArray[i].getPosition(), e, b);
        pusher.push(particleArray[i], e, b, dt);
    }
}

// Reflective boundary conditions
template<class Ensemble>
void applyBoundaryConditions(Ensemble& particles, int beginIdx, int endIdx)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle> particleArray = particles.getParticles();
    pica::Vector3<double> minPosition = particles.getMinPosition();
    pica::Vector3<double> maxPosition = particles.getMaxPosition();
    for (int i = beginIdx; i < endIdx; i++) {
        pica::Vector3<double> position = particleArray[i].getPosition();
        pica::Vector3<double> momentum = particleArray[i].getMomentum();
        for (int d = 0; d < 3; d++)
            if (position[d] < minPosition[d]) {
                position[d] = 2.0 * minPosition[d] - position[d];
                momentum[d] = -momentum[d];
            }
            else
                if (position[d] > maxPosition[d]) {
                    position[d] = 2.0 * maxPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
        particleArray[i].setPosition(position);
        particleArray[i].setMomentum(momentum);
    }
}

template<class Grid>
void updateField(Grid& fields, double dt)
{
    pica::YeeSolver solver;
    solver.updateB(fields, dt / 2.0);
    solver.updateE(fields, dt);
    solver.updateB(fields, dt / 2.0);
}
