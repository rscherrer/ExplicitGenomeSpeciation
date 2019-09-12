#include "MetaPop.h"

double Xst(const vecDbl &v, const vecUns &n)
{
    return 1.0 - (n[0u] * v[0u] - n[1u] * v[1u]) / (n[2u] * v[2u]);
}

size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    StreamBag out;

    if (record)
        out.openAll();

    for (; t < tmax; ++t) {

        // Sort out the sexes
        pops[0u].sortSexes();
        pops[1u].sortSexes();


        // Analyze and record
        if (record && t % tsave == 0u) {

            // Prepare analysis for output

            // Analysis step

            // Classify ecotypes
            // Calculate global statistics
            // Calculate ecotype statistics
            // Calculate divergence statistics
            // Write to output

            //classifyEcotypes();


            // Assign individuals to a group based on trait value

            // Reset statistics
            meanPhenotypes = { zeros(3u), zeros(3u), zeros(3u) };
            meanGenValues = zeros(3u);
            pheVariances = { zeros(3u), zeros(3u), zeros(3u) };
            genVariances = zeros(3u);
            addVariances = zeros(3u);
            domVariances = zeros(3u);
            intVariances = zeros(3u);
            Pst = zeros(3u);

            // Mean phenotypes at the scale of the metapopulation
            size_t metapopsize = 0u;
            for (size_t p = 0u; p < 2u; ++p) {
                for (auto ind : pops[p].individuals) {

                    const vecDbl traitValues = ind->getTraits();
                    const vecDbl geneticValues = ind->getGeneticValues();
                    for (size_t trait = 0u; trait < 3u; ++trait) {
                        meanPhenotypes[trait][2u] += traitValues[trait];
                        meanGenValues[trait] += geneticValues[trait];
                        pheVariances[trait][2u] += sqr(traitValues[trait]);
                        genVariances[trait] += sqr(geneticValues[trait]);
                    }
                }
                metapopsize += pops[p].getPopSize();
            }
            for (size_t trait = 0u; trait < 2u; ++trait) {
                meanPhenotypes[trait][2u] /= metapopsize;
                pheVariances[trait][2u] /= metapopsize;
                pheVariances[trait][2u] -= sqr(meanPhenotypes[trait][2u]);
                genVariances[trait] /= metapopsize;
                genVariances[trait] -= sqr(meanGenValues[trait]);
            }

            // Assign ecotypes and calculate ecotype-specific means
            for (size_t eco = 0u; eco < 2u; ++eco)
                ecotypes[eco].clear();
            for (size_t p = 0u; p < 2u; ++p) {
                for (auto ind : pops[p].individuals) {
                    size_t group = ind->getEcoTrait() < meanPhenotypes[0u][2u];
                    vecDbl traitValues = ind->getTraits();
                    for (size_t trait = 0u; trait < 3u; ++trait)
                        meanPhenotypes[trait][group] += traitValues[trait];
                    ecotypes[group].push_back(ind);
                }
            }

            for (size_t trait = 0u; trait < 2u; ++trait)
                for (size_t eco = 0u; eco < 2u; ++eco)
                    meanPhenotypes[trait][eco] /= ecotypes[eco].size();

            // Locus-specific variance decomposition
            for (size_t locus = 0u; locus < genome.nloci; ++locus) {

                double meanAlleleCount = 0.0;
                double varAlleleCount = 0.0;
                double meanLocusGenValue = 0.0;
                double covGenValueAlleleCount = 0.0;

                vecUns genotypeCounts = uzeros(3u);
                vecDbl meanGenotypeGenValues = zeros(3u);

                for (size_t eco = 0u; eco < 2u; ++eco) {
                    for (auto ind : ecotypes[eco]) {

                        size_t zyg = ind->getZygosity(locus);
                        double genvalue = ind->getLocusGenValue(locus);

                        meanAlleleCount += zyg;
                        varAlleleCount += sqr(zyg);

                        ++genotypeCounts[zyg];
                        meanGenotypeGenValues[zyg] += genvalue;

                        meanLocusGenValue += genvalue;
                        covGenValueAlleleCount += zyg * genvalue;


                    }
                }

                meanAlleleCount /= metapopsize;
                varAlleleCount /= metapopsize;
                varAlleleCount -= meanAlleleCount;
                meanLocusGenValue /= metapopsize;
                covGenValueAlleleCount /= metapopsize;
                covGenValueAlleleCount -= meanAlleleCount * meanLocusGenValue;
                for (size_t zyg = 0u; zyg < 3u; ++zyg)
                    meanGenotypeGenValues[zyg] /= genotypeCounts[zyg];

                const size_t trait = genome.traits[locus];

                double avgMutEffect = covGenValueAlleleCount / varAlleleCount;
                double locusVarA = sqr(avgMutEffect) * varAlleleCount;
                addVariances[trait] += locusVarA;

                vecDbl domDeviations = zeros(3u);
                double locusVarD = 0.0;
                for (size_t zyg = 0u; zyg < 3u; ++zyg) {
                    double breedingValue = avgMutEffect;
                    breedingValue *= (zyg - meanAlleleCount);
                    double addExpectation = meanLocusGenValue - breedingValue;
                    domDeviations[zyg] = meanGenotypeGenValues[zyg];
                    domDeviations[zyg] -= addExpectation;
                    locusVarD += genotypeCounts[zyg] * sqr(domDeviations[zyg]);
                }
                locusVarD /= metapopsize;
                domVariances[trait] += locusVarD;

                double locusVarI = 0.0;
                for (size_t eco = 0u; eco < 2u; ++eco) {
                    for (auto ind : ecotypes[eco]) {
                        double intDeviation = ind->getLocusGenValue(locus);
                        const size_t zyg = ind->getZygosity(locus);
                        intDeviation -= meanGenotypeGenValues[zyg];
                        locusVarI += sqr(intDeviation);
                    }
                }
                locusVarI /= metapopsize;
                locusVarI *= 2.0;

                intVariances[trait] += locusVarI;

            }

            const size_t n0 = ecotypes[0u].size();
            const size_t n1 = ecotypes[1u].size();
            const vecUns census = { n0, n1, metapopsize };
            for (size_t trait = 0u; trait < 3u; ++trait) {
                Pst[trait] = Xst(pheVariances[trait], census);
            }

            // Load output to buffer
            loadBuffer(t);

            // Write to files
            // for (size_t f = 0u; f < out.names.size(); ++f)
                // buffer.write(out.files[f], buffer.fields[f]);

        }

        // Dispersal
        Crowd migrants1 = pops[0u].emigrate(dispersal);
        Crowd migrants2 = pops[1u].emigrate(dispersal);
        pops[0u].immigrate(migrants2);
        pops[1u].immigrate(migrants1);

        // Resource acquisition
        pops[0u].consume();
        pops[1u].consume();

        // Reproduction
        pops[0u].reproduce(birth, sexsel, genome, networks);
        pops[1u].reproduce(birth, sexsel, genome, networks);

        // Survival
        if (!pops[0u].survive(survival) && !pops[1u].survive(survival)) {
            std::cout << "The population went extinct at t = " << t << '\n';
            break;
        }
    }

    if (record)
        out.closeAll();

    return t;
}

void MetaPop::loadBuffer(const size_t &t)
{
    buffer.flush();
    buffer.add(size2dbl(t));

    // Census
    buffer.add(size2dbl(ecotypes[0u].size()));
    buffer.add(size2dbl(ecotypes[1u].size()));
    buffer.add(size2dbl(pops[0u].getPopSize()));
    buffer.add(size2dbl(pops[1u].getPopSize()));
    buffer.add(size2dbl(pops[0u].getNFemales()));
    buffer.add(size2dbl(pops[1u].getNFemales()));

    // Resources in each habitat
    buffer.add(pops[0u].getResources()[0u]);
    buffer.add(pops[1u].getResources()[0u]);
    buffer.add(pops[0u].getResources()[1u]);
    buffer.add(pops[1u].getResources()[1u]);

    // Phenotypes
    for (size_t trait = 0u; trait < 2u; ++trait) {
        for (size_t group = 0u; group < 3u; ++group)
            buffer.add(meanPhenotypes[trait][group]);
        buffer.add(pheVariances[trait][2u]);
        buffer.add(genVariances[trait]);
        buffer.add(addVariances[trait]);
        buffer.add(domVariances[trait]);
        buffer.add(intVariances[trait]);
        buffer.add(Pst[trait]);
    }

    // buffer.add(getEcoIsolation());
    // buffer.add(getSpatialIsolation());
    // buffer.add(getMatingIsolation());
}

double MetaPop::getEcoIsolation(const double &ecomean)
{
    // Ecological isolation is the standard deviation in ecological trait value

    double ei = 0.0;

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t i = 0u; i < pop.getPopSize(); ++i) {
            auto ind = pop.individuals[i];
            ei += sqr(ind->getEcoTrait());
        }
    }

    const size_t n0 = pops[0u].getPopSize();
    const size_t n1 = pops[1u].getPopSize();

    ei /= (n0 + n1);
    ei -= sqr(ecomean); // minus the square of the mean
    return sqrt(ei); // return standard deviation
}

double MetaPop::getSpatialIsolation()
{

    // Calculated as some non-overlap between ecotypes

    // Ecotype-by-habitat table
    std::vector<vecUns> n = { uzeros(2u), uzeros(2u) };

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t i = 0u; i < pop.getPopSize(); ++i) {
            auto ind = pop.individuals[i];
            size_t ecotype = ind->getEcotype();
            ++n[p][ecotype];
        }
    }

    double si = n[0u][0u] * n[1u][1u] - n[0u][1u] * n[1u][0u];
    si /= sqrt(n[0u][0u] + n[0u][1u]);
    si /= sqrt(n[1u][0u] + n[1u][1u]);
    si /= sqrt(n[0u][0u] + n[1u][0u]);
    si /= sqrt(n[0u][1u] + n[1u][1u]);

    return si;
}


double MetaPop::getMatingIsolation()
{

    // Count homogamic and heterogamic crossings

    // For each of them sample a number of males to encounter, from the metapop
    // Evaluate each male by a yes or no
    // Update the table of matings accordingly by looking at ecotypes
    // RI = 0 is mating tests are not possible

    size_t nfemales = 0u;
    size_t nmales = 0u;

    for (size_t p = 0u; p < 2u; ++p) {
        nfemales += pops[p].getNFemales();
        nmales += pops[p].getNMales();
    }

    if (nfemales == 0u || nmales == 0u)
        return 0.0;

    // Table of crossings
    std::vector<vecUns> m = { uzeros(2u), uzeros(2u) };

    // Make a vector with all males of the metapop
    Crowd allMales;
    for (size_t p = 0u; p < 2u; ++p)
        for (size_t i = 0u; i < pops[p].getNMales(); ++i)
            allMales.push_back(pops[p].individuals[i]);

    // Loop through the females of the metapop and test their preference
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t i = 0u; i < pops[p].getNFemales(); ++i) {
            auto fem = pops[p].females[i];
            size_t nencounters = rnd::poisson(1.0 / matingcost);
            while (nencounters) {
                auto candidate = allMales[rnd::random(allMales.size())];
                if (fem->acceptMate(candidate->getEcoTrait(), sexsel))
                    ++m[fem->getEcotype()][candidate->getEcotype()];
                --nencounters;
            }
        }
    }

    double ri = m[0u][0u] * m[1u][1u] - m[0u][1u] * m[1u][0u];
    ri /= sqrt(m[0u][0u] * m[1u][1u] * m[0u][1u] * m[1u][0u]);
    return ri;

}










