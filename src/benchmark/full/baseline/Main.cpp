#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"

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


void runBenchmark(const utility::FullParameters& parameters);

void runIteration(double dt);


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
    runBenchmark(parameters);
    timer->stop();

    utility::printResult(parameters, timer->getElapsed());

    return 0;
}


// Run the whole benchmark
void runBenchmark(const utility::FullParameters& parameters)
{
    omp_set_num_threads(parameters.numThreads);

    // time step
    const double dt = 0.1 / pica::Constants<double>::c();

    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(dt);
}


// Simulate one time step
void runIteration(double dt)
{
}
