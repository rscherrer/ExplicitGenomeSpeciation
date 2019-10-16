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
        size_t size = 0u;
        for (size_t i = 0u; i < population.size(); ++i) {
            if (population[i].getHabitat() == h) {
                ++size;
            }
        }
        return size;
    }
    double getResource(const size_t &h, const size_t &r) const
    {
        return resources[h][r];
    }
    double getSumFitness() const
    {
        double sum = 0.0;
        for (size_t i = 0u; i < population.size(); ++i) {
            sum += population[i].getFitness();
        }
        return sum;
    }
    double getVarFitness() const
    {
        double sum = 0.0;
        double ssq = 0.0;
        for (size_t i = 0u; i < population.size(); ++i) {
            sum += population[i].getFitness();
            ssq += utl::sqr(population[i].getFitness());
        }
        const size_t n = population.size();
        return ssq / n - utl::sqr(sum / n);
    }
    double getMeanEcoTrait(const size_t &h) const // can be removed
    {
        double mean = 0.0;
        size_t n = 0u;
        for (size_t i = 0u; i < population.size(); ++i) {
            if (population[i].getHabitat() == h) {
                ++n;
                mean += population[i].getEcoTrait();
            }
        }
        mean /= n;
        return mean;
    }
    double getMeanEcotype(const size_t &h) const // can be removed
    {
        double mean = 0.0;
        size_t n = 0u;
        for (size_t i = 0u; i < population.size(); ++i) {
            if (population[i].getHabitat() == h) {
                ++n;
                mean += population[i].getEcotype();
            }
        }
        mean /= n;
        return mean;
    }
    double getMeanMatePref() const // can be removed
    {
        double mean = 0.0;
        for (size_t i = 0u; i < population.size(); ++i) {
            mean += population[i].getMatePref();
        }
        mean /= population.size();
        return mean;
    }
    double getFitness(const size_t &i) const // can be removed
    {
        return population[i].getFitness();
    }
    double getFeeding(const size_t &i, const size_t &r) const // can be removed
    {
        return population[i].getFeeding(r);
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

    Crowd population;
    bool isburnin;

    Matrix resources; // per habitat per resource
    MatUns sexcounts; // per habitat per sex

};

#endif
