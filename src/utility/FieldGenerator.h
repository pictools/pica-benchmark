#ifndef PICA_BENCHMARK_UTILITY_FIELD_GENERATOR_H
#define PICA_BENCHMARK_UTILITY_FIELD_GENERATOR_H


#include "utility/Random.h"

#include "pica/math/Vectors.h"

#include <vector>


namespace utility {


template<class T>
std::vector<pica::Vector3<T> > generateField(int size)
{
    Random random;
    std::vector<pica::Vector3<T> > result(size);
    for (int i = 0; i < result.size(); i++) {
        result[i].x = random.getUniform();
        result[i].y = random.getUniform();
        result[i].z = random.getUniform();
    }
    return result;
}


} // namespace utility


#endif
