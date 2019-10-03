#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Deme.h"
#include "Output.h"
#include "Utilities.h"
#include "Types.h"
#include "Stats.h"
#include <cassert>

typedef std::vector<Deme> vecPop;
typedef std::vector<std::ofstream *> vecStreams;

class MetaPop
{

public:

    MetaPop(const Param &pars, const GenArch &arch, const bool &isburnin) :
        population(populate(pars, arch))
    {

        // There should be the right number of individuals in each habitat
        assert(getDemeSize(0u) == pars.demesizes[0u]);
        assert(getDemeSize(1u) == pars.demesizes[1u]);

        // The population should be of a certain size
        assert(population.size() == pars.getInitPopSize());
    }

    ~MetaPop() {}

    size_t getPopSize(const size_t &p) const { return pops[p].getPopSize(); }
    size_t getDemeSize(const size_t&) const;
    size_t getNFemales(const size_t &p) const { return pops[p].getNFemales(); }
    size_t getNOffspring(const size_t&) const;
    size_t getSumEcotypes(const size_t&) const;
    size_t getSumFemEcotypes() const;
    size_t getEcoCount(const size_t &e) const { return stats.getEcoCount(e); }
    double getResource(const size_t&, const size_t&) const;
    double getEcoIsolation() const { return stats.getEcoIsolation(); }
    double getSpatialIsolation() const { return stats.getSpatialIsolation(); }
    double getMatingIsolation() const { return stats.getMatingIsolation(); }
    double getPst(const size_t &trait) const { return stats.getPst(trait); }
    double getVarP(const size_t&, const size_t&) const;
    double getSsqPhe(const size_t&, const size_t&) const;
    double getSumPhe(const size_t&, const size_t&) const;
    double getSumTrait(const size_t&, const size_t&) const;

    void cycle();

    void resetEcoTraits(const size_t&, const double&);
    void resetMatePrefs(const size_t&, const double&);
    void resetEcotypes(const size_t&, const size_t&);
    void resetGenders(const size_t&, const bool&);

private:

    Crowd populate(const Param&, const GenArch&);

    void disperse();
    void consume();
    void reproduce();
    void survive();

    Crowd population;

};

#endif
