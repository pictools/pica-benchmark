#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

#include "pica/currentDeposition/CurrentDepositor.h"
#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/fieldSolver/YeeSolver.h"
#include "pica/grid/YeeGrid.h"
#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleArray.h"
#include "pica/particlePush/BorisPusher.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>

template<class Ensemble, class Grid>
void runBenchmark(Ensemble& particles, Grid& fields, const utility::FullParameters& parameters);

template<class Ensemble, class Grid>
void runIteration(Ensemble& particles, Grid& fields, std::vector<Grid>& threadFields, double dt);

template<class Ensemble, class Grid>
void updateParticles(Ensemble& particles, const Grid& fields,
    std::vector<Grid>& threadFields, double dt);

template<class Ensemble, class Grid>
void process(Ensemble& particles, const Grid& fields, std::vector<Grid>& threadFields,
    int beginIdx, int endIdx, double dt);

template<class Ensemble>
void applyBoundaryConditions(Ensemble& particles, int beginIdx, int endIdx);

template<class Ensemble, class Grid>
void push(Ensemble& particles, const Grid& fields,
    int beginIdx, int endIdx, double dt);

template<class Grid>
void finalizeCurrents(Grid& fields, std::vector<Grid>& threadFields);

template<class Ensemble, class Grid>
void depositCurrents(Ensemble& particles, std::vector<Grid>& threadFields,
    int beginIdx, int endIdx, double dt);

template<class Grid>
void updateFields(Grid& fields, double dt);

int main(int argc, char* argv[])
{
    utility::FullParameters parameters = utility::readFullParameters(argc, argv);
    utility::printHeader("full-sorted benchmark: using sorted in space 3D particle ensemble, CIC form factor and SoA particle representation",
        parameters);

    // Generate particles randomly,
    // particular coordinates are not important for this benchmark
    typedef pica::Particle<pica::Three> Particle;
    typedef pica::ParticleArraySoA<Particle> ParticleArray;
    typedef typename pica::Ensemble<ParticleArray,
        pica::EnsembleRepresentation_Ordered>::Type Particles;
    Particles particles = utility::generateParticles<Particles>(
        parameters.numCells, parameters.particlesPerCell, parameters.numParticleTypes);

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
    std::vector<Grid> threadFields(parameters.numThreads, fields);

    // time step
    const double dt = 1.0 / (8 * parameters.numCells.x * pica::Constants<double>::c());

    for (int i = 0; i < parameters.numIterations; i++) {
        if (i % parameters.sortingPeriod)
            particles.reorder();
        runIteration(particles, fields, threadFields, dt);
    }
}


// Simulate one time step
template<class Ensemble, class Grid>
void runIteration(Ensemble& particles, Grid& fields,
    std::vector<Grid>& threadFields, double dt)
{
    updateParticles(particles, fields, threadFields, dt);
    finalizeCurrents(fields, threadFields);
    updateFields(fields, dt);
}

template<class Ensemble, class Grid>
void updateParticles(Ensemble& particles, const Grid& fields,
    std::vector<Grid>& threadFields, double dt)
{
    // Each thread processes some particles
    const int numParticles = particles.size();
    const int numThreads = pica::getNumThreads();
    const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
    #pragma omp parallel for
    for (int idx = 0; idx < numThreads; idx++) {
        const int beginIdx = idx * particlesPerThread;
        const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
        process(particles, fields, threadFields, beginIdx, endIdx, dt);
    }
}

template<class Ensemble, class Grid>
void process(Ensemble& particles, const Grid& fields, std::vector<Grid>& threadFields,
    int beginIdx, int endIdx, double dt)
{
    push(particles, fields, beginIdx, endIdx, dt);
    applyBoundaryConditions(particles, beginIdx, endIdx);
    depositCurrents(particles, threadFields, beginIdx, endIdx, dt);
}

template<class Ensemble, class Grid>
void push(Ensemble& particles, const Grid& fields,
    int beginIdx, int endIdx, double dt)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles();
    pica::BorisPusher<Particle> pusher;
    pica::FieldInterpolatorCIC<Grid> fieldInterpolator(fields);
    #pragma omp simd
    #pragma forceinline
    for (int i = beginIdx; i < endIdx; i++) {
        pica::Vector3<double> e, b;
        fieldInterpolator.get(particleArray[i].getPosition(), e, b);
        pusher.push(&particleArray[i], e, b, dt);
    }
}

// Reflective boundary conditions
template<class Ensemble>
void applyBoundaryConditions(Ensemble& particles, int beginIdx, int endIdx)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles();
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

template<class Ensemble, class Grid>
void depositCurrents(Ensemble& particles, std::vector<Grid>& threadFields,
    int beginIdx, int endIdx, double dt)
{
    Grid& fields = threadFields[omp_get_thread_num()];
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles();
    const double halfDt = 0.5 * dt;
    pica::CurrentDepositorCIC<Grid> currentDepositor(fields);

    // Zeroise currents
    for (int i = 0; i < fields.getSize().x; i++)
    for (int j = 0; j < fields.getSize().y; j++)
    for (int k = 0; k < fields.getSize().z; k++) {
        fields.jx(i, j, k) = 0.0;
        fields.jy(i, j, k) = 0.0;
        fields.jz(i, j, k) = 0.0;
    }

    for (int i = beginIdx; i < endIdx; i++) {
        pica::Vector3<double> position = particleArray[i].getPosition();
        for (int d = 0; d < 3; d++)
            position[d] -= particleArray[i].getVelocity()[d] * halfDt;
        pica::Vector3<double> current = particleArray[i].getVelocity() *
            particleArray[i].getCharge() * (double)particleArray[i].getFactor();
        currentDepositor.deposit(position, current);
    }
}

template<class Grid>
void finalizeCurrents(Grid& fields, std::vector<Grid>& threadFields)
{
    double normalization = 1.0 / fields.getStep().volume();
    const int sizeX = fields.getSize().x;
    const int sizeY = fields.getSize().y;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < sizeX; i++)
    for (int j = 0; j < sizeY; j++)
    for (int k = 0; k < fields.getSize().z; k++) {
        fields.jx(i, j, k) = 0.0;
        fields.jy(i, j, k) = 0.0;
        fields.jz(i, j, k) = 0.0;
        for (int threadIdx = 0; threadIdx < threadFields.size(); threadIdx++) {
            const Grid& currentGrid = threadFields[threadIdx];
            fields.jx(i, j, k) += currentGrid.jx(i, j, k);
            fields.jy(i, j, k) += currentGrid.jy(i, j, k);
            fields.jz(i, j, k) += currentGrid.jz(i, j, k);
        }
        fields.jx(i, j, k) *= normalization;
        fields.jy(i, j, k) *= normalization;
        fields.jz(i, j, k) *= normalization;
    }
}

template<class Grid>
void updateFields(Grid& fields, double dt)
{
    pica::YeeSolver solver;
    solver.updateB(fields, dt / 2.0);
    solver.updateE(fields, dt);
    solver.updateB(fields, dt / 2.0);
}
