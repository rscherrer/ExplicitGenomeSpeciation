#include "Pedigree.h"

// Constructor
Pedigree::Pedigree(const std::string &pedigreename, const bool &pedigreesave) :
    pedigree(new std::ofstream)
{

    if (pedigreesave) {

        // Open the pedigree file
        pedigree->open(pedigreename, std::ios::binary);
        if (!pedigree->is_open()) {
            std::string msg = "Unable to open output pedigree file";
            throw std::runtime_error(msg);
        }

    }

}

// Destructor
Pedigree::~Pedigree()
{
    shutdown();
}

void Pedigree::shutdown() {
    pedigree->close();
}

// Member functions
//-----------------

void Pedigree::analyze(const MetaPop &m, const Param &p, const GenArch &arch) {

    // Analyze the distributions of offspring

    // Sort males and females in the population
    std::vector<size_t> males;
    std::vector<size_t> females;
    males.reserve(m.population.size());
    females.reserve(m.population.size());
    for (size_t i =0u; i < m.population.size(); ++i) {
        const size_t sex = m.population[i].getGender();
        if (sex) females.push_back(i); else males.push_back(i);
    }
    males.shrink_to_fit();
    females.shrink_to_fit();

    if (females.size() && males.size()) {

        size_t ntrials = p.pedigreetrials;

        // Sample from a distribution of males and a distribution of females
        auto femalepool = rnd::random(0u, females.size() - 1u);
        auto malepool = rnd::random(0u, males.size() - 1u);

        // Sample many pairs of males and females with replacement
        while (ntrials) {

            const size_t fem = females[femalepool(rnd::rng)];
            const size_t mal = males[malepool(rnd::rng)];
            const size_t ecof = m.population[fem].getEcotype();
            const size_t ecom = m.population[mal].getEcotype();

            stf::write(static_cast<double>(fem), pedigree);
            stf::write(static_cast<double>(mal), pedigree);
            stf::write(static_cast<double>(ecof), pedigree);
            stf::write(static_cast<double>(ecom), pedigree);

            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(m.population[fem].getTraitValue(trait), pedigree);

            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(m.population[mal].getTraitValue(trait), pedigree);

            // Produce a fixed number of offspring
            for (size_t k = 0u; k < p.pedigreeoffspring; ++k) {

                const Individual baby = Individual(p, arch, m.population[fem], m.population[mal]);

                for (size_t trait = 0u; trait < 3u; ++trait)
                    stf::write(baby.getTraitValue(trait), pedigree);

            }

            --ntrials;
        }
    }
}
