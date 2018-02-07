#ifndef PICA_BENCHMARK_UTILITY_RANDOM_H
#define PICA_BENCHMARK_UTILITY_RANDOM_H


namespace utility {


// This class implements the same random number generator that the NAG
// library does.
// See http://nag.co.uk/numeric/CL/nagdoc_cl09/xhtml/G05/g05cac.xml
class Random
{
    static const unsigned long long modulo = 1ULL << 59; /* 2^59 */
    static const unsigned long long mask = modulo - 1; /* for taking the remainder */
    static const unsigned long long default_multiplier = 302875106592253ULL; /* 13^13 */

    unsigned long long multiplier;
    unsigned long long state;

public:

    Random(unsigned long long seed = 0);

    /* Get a normally distributed random number with zero mean and
       a variance of 1. This requires two elements from the random stream. */
    double getNormal();
    /* Get a uniformly distributed random number in [0.0; 1.0). */
    double getUniform();
    /* Consume every nth element of the random stream. */
    void setLeapfrog(int n);
    /* Skip n elements of the random stream. */
    void skipAhead(int n);
};


} // namespace utility


#endif
