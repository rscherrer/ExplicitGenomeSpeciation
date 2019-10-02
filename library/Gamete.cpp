#include "Gamete.h"

void Gamete::mutate(const double &rate)
{

    // Sample a number of mutations from a poisson
    // Sample the mutated targets
    // Flip the alleles

    size_t nmut = rnd::poisson(rate * seq.size());

    while (nmut) {

        size_t target = rnd::random(seq.size());
        seq[target] = seq[target] ? 0u : 1u;

        --nmut;
    }

}


