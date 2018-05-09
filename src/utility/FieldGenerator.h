#ifndef PICA_BENCHMARK_UTILITY_FIELD_GENERATOR_H
#define PICA_BENCHMARK_UTILITY_FIELD_GENERATOR_H


#include "utility/Random.h"

#include "pica/math/Vectors.h"

#include <vector>


namespace utility {


template<class T>
pica::Vector3<T> generateField()
{
    Random random;
    pica::Vector3<T> result;
    result.x = random.getUniform();
    result.y = random.getUniform();
    result.z = random.getUniform();
    return result;
}


template<class Grid>
Grid generateField(typename Grid::PositionType minPosition,
    typename Grid::PositionType maxPosition,
    typename Grid::IndexType numInternalCells)
{
    typedef typename Grid::PositionType PositionType;
    typedef typename Grid::IndexType IndexType;
    PositionType step = (maxPosition - minPosition) / (typename Grid::PositionType(numInternalCells));
    int numGhostCells = 2;
    PositionType origin = minPosition - step * static_cast<typename pica::ScalarType<PositionType>::Type>(numGhostCells);
    IndexType numCells = numInternalCells;
    for (int d = 0; d < pica::VectorDimensionHelper<IndexType>::dimension; d++)
        numCells[d] += 2 * numGhostCells;
    Grid fields(origin, step, numCells);

    Random random;
    IndexType size = fields.getSize();
    for (int i = 0; i < size.x; i++)
    for (int j = 0; j < size.y; j++)
    for (int k = 0; k < size.z; k++) {
        fields.ex(i, j, k) = random.getUniform();
        fields.ey(i, j, k) = random.getUniform();
        fields.ez(i, j, k) = random.getUniform();
        fields.bx(i, j, k) = random.getUniform();
        fields.by(i, j, k) = random.getUniform();
        fields.bz(i, j, k) = random.getUniform();
    }
    return fields;
}


} // namespace utility


#endif
