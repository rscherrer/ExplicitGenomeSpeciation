#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Param.h"
#include "GenArch.h"
#include "Individual.h"
#include "Utilities.h"
#include "Types.h"
#include <cassert>

typedef std::vector<Individual> Crowd;

class MetaPop
{

    friend class Collector;

public:

    MetaPop(const Param &pars, const GenArch &arch) :
        population(populate(pars, arch)),
        isburnin(pars.tburnin > 0),
        resources(utl::zeros(2u, 2u)),
        sexcounts(utl::uzeros(2u, 2u))
    {

        // Test right number of individuals
        // Test right number of individuals in each habitat

    }

    ~MetaPop() {}

    void cycle(const Param&, const GenArch&);
    void exitburnin();
    bool isextinct() const;

    // Getters called from outside
    size_t getSize() const
    {
        return population.size();
    }
    size_t getDemeSize(const size_t &h) const
    {
        return sexcounts[h][0u] + sexcounts[h][1u];
    }
    double getResource(const size_t &h, const size_t &r) const
    {
        return resources[h][r];
    }

    // Resetters used in tests
    void resetEcoTraits(const double &x, const Param &p)
    {
        for (size_t i = 0u; i < population.size(); ++i){
            population[i].resetEcoTrait(x, p);
        }
    }
    void resetEcoTraits(const size_t &h, const double &x, const Param &p)
    {
        for (size_t i = 0u; i < population.size(); ++i) {
            if (population[i].getHabitat() == h) {
                population[i].resetEcoTrait(x, p);
            }
        }
    }
    void resetMatePrefs(const double &x)
    {
        for (size_t i = 0u; i < population.size(); ++i) {
            population[i].resetMatePref(x);
        }
    }
    void resetGenders(const bool &sex)
    {
        for (size_t i = 0u; i < population.size(); ++i) {
            population[i].resetGender(sex);
        }
    }

private:

    Crowd populate(const Param&, const GenArch&);

    void disperse(const Param&);
    void consume(const Param&);
    void reproduce(const Param&, const GenArch&);
    void survive(const Param&);

    MatUns matingtrials(const Param&) const;

    Crowd population;
    bool isburnin;

    Matrix resources; // per habitat per resource
    MatUns sexcounts; // per habitat per sex

};

#endif
