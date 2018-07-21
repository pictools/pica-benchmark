#ifndef PICA_PARTICLE_BASELINE_H
#define PICA_PARTICLE_BASELINE_H


#include "pica/math/Constants.h"
#include "pica/math/FP.h"
#include "pica/math/Vectors.h"

#include <cmath>
#include <vector>


namespace pica {

typedef FP Real;
typedef Real MassType;
typedef Real ChargeType;

struct ParticleType {
    MassType mass;
    ChargeType charge;
};

namespace ParticleTypes
{
    extern std::vector<ParticleType> typesVector;
    extern const ParticleType* types;
    extern int numTypes;
};

// Baseline representation for particle
template<Dimension dimension>
class ParticleBaseline {
public:

    // Types for conforming ParticleInterface
    typedef FP Real;
    typedef typename VectorTypeHelper<dimension, Real>::Type PositionType;
    typedef typename VectorTypeHelper<Three, Real>::Type MomentumType;
    typedef Real GammaType;
    typedef Real MassType;
    typedef Real ChargeType;
    typedef Real FactorType;
    typedef short TypeIndexType;

    ParticleBaseline() :
        factor(1),
        typeIndex(0)
    {}

    ParticleBaseline(const PositionType& position, const MomentumType& momentum, FactorType factor = 1, TypeIndexType typeIndex = 0) :
        position(position), momentum(momentum), mass(mass), charge(charge), factor(factor), typeIndex(typeIndex) {}

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const
    {
        return p * Constants<MassType>::c() * getMass();
    }

    void setMomentum(const MomentumType& newMomentum)
    {
        p = newMomentum / (Constants<GammaType>::c() * getMass());
    }

    MomentumType getP() const { return p; }
    void setP(const MomentumType& newP) { p = newP; }

    MomentumType getVelocity() const
    {
        return p * Constants<GammaType>::c() / sqrt((FP)1 + p.norm2());
    }

    void setVelocity(const MomentumType& newVelocity)
    {
        p = newVelocity / sqrt(constants::c * constants::c - newVelocity.norm2());
    }

    GammaType getGamma() const { return sqrt(static_cast<FP>(1.0) + p.norm2()); }

    MassType getMass() const { return ParticleTypes::types[typeIndex].mass; }

    ChargeType getCharge() const { return ParticleTypes::types[typeIndex].charge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

    TypeIndexType getType() const { return typeIndex; }
    void setType(TypeIndexType newType) { typeIndex = newType; }

private:

    PositionType position;
    MomentumType p;
    FactorType factor;
    TypeIndexType typeIndex;
};


} // namespace pica


#endif