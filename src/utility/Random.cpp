#include "utility/Random.h"

#include "pica/math/Constants.h"

#include <cmath>
#include <limits>


namespace {

bool isNumberFinite(double x)
{
    return x <= std::numeric_limits<double>::max()
        && x >= -std::numeric_limits<double>::max();
}

} // anonymous namespace

namespace utility {


Random::Random(unsigned long long seed)
{
    state = seed * 2 + 1;
    multiplier = default_multiplier;
}


double Random::getNormal()
{
    /* This is the standard Box-Muller transform, except we only use
       one of the results, for simplicity. */
    double result;
    do
    {
        double u1 = getUniform();
        double u2 = getUniform();
        result = std::sqrt(-2 * std::log(u1)) * std::cos(2 * pica::Constants<double>::pi() * u2);
    } while (!isNumberFinite(result));
    return result;
}


double Random::getUniform()
{
    state *= multiplier;
    return (state & mask) / double(modulo);
}


unsigned long long ullpow(unsigned long long base,
                          unsigned long long exponent) {
    unsigned long long result = 1;
    while (exponent != 0) {
        if (exponent & 1)
            result *= base;
        base *= base;
        exponent >>= 1;
    }
    return result;
}


void Random::setLeapfrog(int n)
{
    multiplier = ullpow(default_multiplier, n);
}


void Random::skipAhead(int n)
{
    state *= ullpow(multiplier, n);
}


} // namespace utility
