#include "utility/FieldGenerator.h"
#include "utility/Output.h"
#include "utility/Parameters.h"
#include "utility/Parser.h"
#include "utility/ParticleGenerator.h"
#include "utility/Random.h"
#include "utility/Timer.h"
#include "pica/math/Constants.h"
#include "pica/particles/Particle.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>

using namespace pica;

int main(int argc, char* argv[])
{
    utility::PusherParameters parameters = utility::readPusherParameters(argc, argv);
    utility::printHeader("pusher-vectorized benchmark: using optimized 3D Boris particle pusher implementation with precomputing of inverse gamma and SoA particle representation",
        parameters);

    std::vector<double> dataX(parameters.numParticles);
    std::vector<double> dataY(parameters.numParticles);
    std::vector<double> dataZ(parameters.numParticles);
    std::vector<double> dataPx(parameters.numParticles);
    std::vector<double> dataPy(parameters.numParticles);
    std::vector<double> dataPz(parameters.numParticles);
    std::vector<double> dataCoeff(parameters.numParticleTypes);
    std::vector<int> dataTypeIndex(parameters.numParticles);
    double* x = &dataX[0];
    double* y = &dataY[0];
    double* z = &dataZ[0];
    double* px = &dataPx[0];
    double* py = &dataPy[0];
    double* pz = &dataPz[0];
    double* coeff = &dataCoeff[0];
    int* typeIndex = &dataTypeIndex[0];

    // Generate particles randomly,
    // particular coordinates and other data are not important for this benchmark
    utility::Random random;
    utility::detail::initParticleTypes(parameters.numParticleTypes);
    for (int i = 0; i < parameters.numParticles; i++)
    {
        pica::Particle3d particle;
        utility::detail::generateParticle(particle, random, parameters.numParticleTypes);

        x[i] = particle.getPosition().x;
        y[i] = particle.getPosition().y;
        z[i] = particle.getPosition().z;
        px[i] = particle.getP().x;
        py[i] = particle.getP().y;
        pz[i] = particle.getP().z;
        typeIndex[i] = particle.getType();
    }

    typedef pica::Vector3<double> FieldValue;
    FieldValue electricFieldValue = utility::generateField<double>();
    FieldValue magneticFieldValue = utility::generateField<double>();

    double dt = (rand() / (double)RAND_MAX) * 1e-10;
    for (int i = 0; i < ParticleTypes::numTypes; i++)
        coeff[i] = ParticleTypes::types[i].charge * dt / (2.0 * ParticleTypes::types[i].mass * pica::Constants<double>::c());

    std::auto_ptr<utility::Stopwatch> timer(utility::createStopwatch());

    timer->start();
    for (int iter = 0; iter < parameters.numIterations; iter++)
    {
        const int numThreads = pica::getNumThreads();
        const int chunkSize = (parameters.numParticles + numThreads - 1) / numThreads;

        #pragma omp parallel for
        for (int idx = 0; idx < numThreads; idx++)
        {
            const int beginIdx = idx * chunkSize;
            const int endIdx = std::min(beginIdx + chunkSize, parameters.numParticles);

            #pragma omp simd
            for (int i = beginIdx; i < endIdx; i++)
            {
                double eCoeff = coeff[typeIndex[idx]];
                double eMomentumX = electricFieldValue.x * eCoeff;
                double eMomentumY = electricFieldValue.y * eCoeff;
                double eMomentumZ = electricFieldValue.z * eCoeff;
                double umX = px[i] + eMomentumX;
                double umY = py[i] + eMomentumY;
                double umZ = pz[i] + eMomentumZ;
                double cf = eCoeff / sqrt(1.0 + umX * umX + umY * umY + umZ * umZ);
                double tX = magneticFieldValue.x * cf;
                double tY = magneticFieldValue.y * cf;
                double tZ = magneticFieldValue.z * cf;
                double uprimeX = umX + umZ * tZ - umZ * tY;
                double uprimeY = umY + umY * tX - umX * tZ;
                double uprimeZ = umZ + umX * tY - umY * tX;
                cf = 2.0 / (1.0 + tX * tX + tY * tY + tZ * tZ);
                double sX = tX * cf;
                double sY = tY * cf;
                double sZ = tZ * cf;
                px[i] = umX + uprimeY * sZ - uprimeZ * sY + eMomentumX;
                py[i] = umY + uprimeZ * sX - uprimeX * sZ + eMomentumY;
                pz[i] = umZ + uprimeX * sY - uprimeY * sX + eMomentumZ;
                double invGamma = 1.0 / sqrt(1.0 + px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);

                cf = pica::Constants<double>::c() * invGamma * dt;
                x[i] += px[i] * cf;
                y[i] += py[i] * cf;
                z[i] += pz[i] * cf;
            }
        }
    }

    timer->stop();
    utility::printResult(parameters, timer->getElapsed());

    return 0;
}
