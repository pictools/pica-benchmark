#ifndef PICA_BENCHMARK_UTILITY_PARAMETERS_H
#define PICA_BENCHMARK_UTILITY_PARAMETERS_H


#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleArray.h"


namespace utility {


// Parameters for particle pusher benchmarks
struct PusherParameters {
    int numParticles;
    int numParticleTypes;
    int numIterations;
    int numThreads;
};

// Parameters for full particle-in-cell benchmarks
struct FullParameters {
    pica::Vector3<int> numCells;
    int numIterations;
    int particlesPerCell;
    int numParticleTypes;
    double temperature;
    pica::ParticleRepresentation particleRepresentation;
    pica::EnsembleRepresentation ensembleRepresentation;
    int sortingPeriod; // used only for ordered ensemble representation
    pica::Vector3<int> numCellsPerSupercell; // used only for supercell ensemble representation
    bool enablePreloading; // used only for supercell ensemble representation
    int tileSize;
    int numThreads;
};


} // namespace utility


#endif
