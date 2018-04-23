#include "utility/Output.h"

#include "utility/Parameters.h"

#include "pica/threading/OpenMPHelper.h"

#include <iostream>

using namespace std;
using namespace pica;


namespace utility {


void printHeader(const string& message, const FullParameters& parameters)
{
    if (!message.empty())
        cout << message << "\n";
    cout << "Parameters:\n";
    cout << "    Grid size: ";
    for (int d = 0; d < 2; d++)
        cout << parameters.numCells[d] << "x";
    cout << parameters.numCells[2] << "\n";
    cout << "    Number of time iterations: " << parameters.numIterations << "\n";
    cout << "    Particles per cell: " << parameters.particlesPerCell << "\n";
    cout << "    Particles temperature: " << parameters.temperature << "\n";
    cout << "    Particle representation: " << toString(parameters.particleRepresentation) << "\n";
    cout << "    Ensemble representation: " << toString(parameters.ensembleRepresentation);
    if (parameters.ensembleRepresentation == EnsembleRepresentation_Ordered)
        cout << ", sorting period = " << parameters.sortingPeriod;
    else if (parameters.ensembleRepresentation == EnsembleRepresentation_Supercells)
    {
        cout << ", supercell size = " << parameters.numCellsPerSupercell.x << "x" <<
            parameters.numCellsPerSupercell.y << "x" << parameters.numCellsPerSupercell.z;
        cout << ", preloading ";
        if (parameters.enablePreloading)
            cout << "enabled";
        else
            cout << "disabled";
    }
    cout << "\n";
    if (parameters.tileSize)
        cout << "    Tile size: " << parameters.tileSize << "\n";
    if (pica::useOpenMP())
        cout << "    " << parameters.numThreads << " OpenMP threads\n";
    else
        cout << "    OpenMP is disabled\n";
    cout << endl;
}

void printHeader(const string& message, const PusherParameters& parameters)
{
    if (!message.empty())
        cout << message << "\n";
    cout << "Parameters:\n";
    cout << "    Number of particles: " << parameters.numParticles << "\n";
    cout << "    Number of time iterations: " << parameters.numIterations << "\n";
    if (pica::useOpenMP())
        cout << "    " << parameters.numThreads << " OpenMP threads\n";
    else
        cout << "    OpenMP is disabled\n";
    cout << endl;
}


void printResult(const FullParameters& parameters, double runTimeSeconds)
{
    cout << "Total run time: " << runTimeSeconds << " sec";
    long numParticles = parameters.numCells.volume() * parameters.particlesPerCell;
    long numParticleUpdates = numParticles * parameters.numIterations;
    double nsPerParticleUpdate = runTimeSeconds * 1e9 / numParticleUpdates;
    cout << ", " << nsPerParticleUpdate << " ns per particle update";
    cout << endl;
}

void printResult(const PusherParameters& parameters, double runTimeSeconds)
{
    cout << "Total run time: " << runTimeSeconds << " sec";
    long numParticleUpdates = parameters.numParticles * parameters.numIterations;
    double nsPerParticleUpdate = runTimeSeconds * 1e9 / numParticleUpdates;
    cout << ", " << nsPerParticleUpdate << " ns per particle update";
    cout << endl;
}


} // namespace utility
