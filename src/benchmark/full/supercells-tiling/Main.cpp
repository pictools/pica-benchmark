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
void runIteration(Ensemble& particles, Grid& fields, Ensemble& migratingParticles, double dt);

template<class Ensemble, class Grid>
void updateParticles(Ensemble& particles, Grid& fields, Ensemble& migratingParticles, double dt);

template<class Ensemble, class Grid>
void process(Ensemble& particles, Grid& fields, Ensemble& migratingParticles,
    pica::Int3 supercellIdx, double dt);

template<class Ensemble>
void migrateAndApplyBoundaryConditions(Ensemble& particles, Ensemble& migratingParticles, pica::Int3 supercellIdx);

template<class Ensemble, class Grid>
void push(Ensemble& particles, const Grid& fields, pica::Int3 supercellIdx, double dt);

template<class Ensemble, class Grid>
void depositCurrents(Ensemble& particles, Grid& fields,
    pica::Int3 supercellIdx, double dt);

template<class Grid>
void updateFields(Grid& fields, double dt);

int main(int argc, char* argv[])
{
    utility::FullParameters parameters = utility::readFullParameters(argc, argv);
    utility::printHeader("full-supercells-tiling benchmark: using supercell-based 3D particle ensemble with tiling, CIC form factor and SoA particle representation",
        parameters);

    // Generate particles randomly,
    // particular coordinates are not important for this benchmark
    typedef pica::Particle<pica::Three> Particle;
    typedef pica::ParticleArraySoA<Particle> ParticleArray;
    typedef typename pica::Ensemble<ParticleArray,
        pica::EnsembleRepresentation_Supercells>::Type Particles;
    Particles particles = utility::generateParticles<Particles>(
        parameters.numCells, parameters.particlesPerCell, parameters.numCells,
        parameters.numCellsPerSupercell, parameters.numParticleTypes);

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
    Ensemble migratingParticles(particles.getMinPosition(), particles.getMaxPosition(),
        particles.getNumCells(), particles.getNumCellsPerSupercell());

    // time step
    const double dt = 1.0 / (8 * parameters.numCells.x * pica::Constants<double>::c());

    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(particles, fields, migratingParticles, dt);
}


// Simulate one time step
template<class Ensemble, class Grid>
void runIteration(Ensemble& particles, Grid& fields, Ensemble& migratingParticles, double dt)
{
    updateParticles(particles, fields, migratingParticles, dt);
    updateFields(fields, dt);
}

template<class Ensemble, class Grid>
void updateParticles(Ensemble& particles, Grid& fields, Ensemble& migratingParticles, double dt)
{
    pica::Int3 numSupercells = particles.getNumSupercells();

    pica::Int3 superCellStep(3, 3, 3);
    for (int d = 0; d < 3; d++)
        if (particles.getNumCellsPerSupercell()[d] == 1)
            superCellStep[d] = 4;
    for (int startI = 0; startI < superCellStep.x; startI++)
        for (int startJ = 0; startJ < superCellStep.y; startJ++)
            for (int startK = 0; startK < superCellStep.z; startK++) {
#pragma omp parallel for collapse(3)
                for (int i = startI; i < numSupercells.x; i += superCellStep.x)
                    for (int j = startJ; j < numSupercells.y; j += superCellStep.y)
                        for (int k = startK; k < numSupercells.z; k += superCellStep.z)
                            process(particles, fields, migratingParticles, pica::Int3(i, j, k), dt);
            }

    // Finalize migration
#pragma omp parallel for collapse(3)
    for (int i = 0; i < numSupercells.x; i++)
        for (int j = 0; j < numSupercells.y; j++)
            for (int k = 0; k < numSupercells.z; k++) {
                pica::ParticleArraySoA<pica::Particle<pica::Three> >& migrating = migratingParticles.getParticles(pica::Int3(i, j, k));
                for (int idx = 0; idx < migrating.size(); idx++)
                    particles.add(migrating[idx]);
                int size = migrating.size();
                for (int idx = 0; idx < size; idx++)
                    migrating.popBack();
            }
}

template<class Ensemble, class Grid>
void process(Ensemble& particles, Grid& fields, Ensemble& migratingParticles,
    pica::Int3 supercellIdx, double dt)
{
    push(particles, fields, supercellIdx, dt);
    migrateAndApplyBoundaryConditions(particles, migratingParticles, supercellIdx);
    depositCurrents(particles, fields, supercellIdx, dt);
}

template<class Ensemble, class Grid>
void push(Ensemble& particles, const Grid& fields, pica::Int3 supercellIdx, double dt)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles(supercellIdx);
    pica::BorisPusher<Particle> pusher;
    pica::FP3 supercellMinPosition = particles.getMinPosition() +
        fields.getStep() * pica::FP3(supercellIdx * particles.getNumCellsPerSupercell());
    pica::FieldInterpolatorCICSupercell<double> fieldInterpolator(fields,
        supercellMinPosition, particles.getNumCellsPerSupercell());
#pragma omp simd
#pragma forceinline
    for (int i = 0; i < particleArray.size(); i++) {
        pica::Vector3<double> e, b;
        fieldInterpolator.get(particleArray[i].getPosition(), e, b);
        pusher.push(&particleArray[i], e, b, dt);
    }
}

// Reflective boundary conditions
template<class Ensemble>
void migrateAndApplyBoundaryConditions(Ensemble& particles, Ensemble& migratingParticles, pica::Int3 supercellIdx)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles(supercellIdx);
    pica::Vector3<double> minPosition = particles.getMinPosition();
    pica::Vector3<double> maxPosition = particles.getMaxPosition();
    for (int i = 0; i < particleArray.size(); i++) {
        // First apply boundary conditions as it could change the position
        // Non-migrating particles are definitely not subject for BC
        if (particles.getSupercellIndex(particleArray[i]) != supercellIdx) {
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
        // Now check again that migration happens after applying boundary conditions
        if (particles.getSupercellIndex(particleArray[i]) != supercellIdx) {
            migratingParticles.add(particleArray[i]);
            typename Ensemble::ParticleRef lastParticle = particleArray.back();
            particleArray[i].setPosition(lastParticle.getPosition());
            particleArray[i].setType(lastParticle.getType());
            particleArray[i].setFactor(lastParticle.getFactor());
            particleArray.popBack();
            i--;
        }
    }
}

template<class Ensemble, class Grid>
void depositCurrents(Ensemble& particles, Grid& fields,
    pica::Int3 supercellIdx, double dt)
{
    typedef typename Ensemble::Particle Particle;
    pica::ParticleArraySoA<Particle>& particleArray = particles.getParticles(supercellIdx);
    const double halfDt = 0.5 * dt;
    pica::FP3 supercellMinPosition = particles.getMinPosition() +
        fields.getStep() * pica::FP3(supercellIdx * particles.getNumCellsPerSupercell());
    pica::CurrentDepositorCICSupercell<double> currentDepositor(fields, supercellMinPosition, particles.getNumCellsPerSupercell());

    for (int i = 0; i < particleArray.size(); i++) {
        pica::Vector3<double> position = particleArray[i].getPosition();
        for (int d = 0; d < 3; d++)
            position[d] -= particleArray[i].getVelocity()[d] * halfDt;
        pica::Vector3<double> current = particleArray[i].getVelocity() *
            particleArray[i].getCharge() * (double)particleArray[i].getFactor();
        currentDepositor.deposit(position, current);
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
