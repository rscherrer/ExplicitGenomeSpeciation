#include "MetaPop.h"

typedef std::vector<vecDbl> Matrix;
typedef std::vector<vecUns> MatUns;

double Xst(const vecDbl &v, const vecUns &n)
{
    return 1.0 - (n[0u] * v[0u] - n[1u] * v[1u]) / (n[2u] * v[2u]);
}

int MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    int t = - tburnin;
    StreamBag out;

    if (record)
        out.openAll();

    for (; t < tmax; ++t) {

        // Sort out the sexes
        pops[0u].sortSexes();
        pops[1u].sortSexes();


        // Analyze and record
        if (record && t % tsave == 0u && t > 0) {

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
            pheVariances = { zeros(3u), zeros(3u), zeros(3u) };
            genVariances = { zeros(3u), zeros(3u), zeros(3u) };
            addVariances = { zeros(3u), zeros(3u), zeros(3u) };
            nadVariances = { zeros(3u), zeros(3u), zeros(3u) };
            domVariances = zeros(3u);
            intVariances = zeros(3u);
            Pst = zeros(3u);
            Gst = zeros(3u);
            Qst = zeros(3u);
            Fst = zeros(3u);
            vecDbl varS = zeros(3u);
            vecDbl varT = zeros(3u);

            PstScan = zeros(genome.nloci);
            GstScan = zeros(genome.nloci);
            QstScan = zeros(genome.nloci);
            CstScan = zeros(genome.nloci);
            FstScan = zeros(genome.nloci);

            Matrix meanGenValues = { zeros(3u), zeros(3u), zeros(3u) };

            // Mean phenotypes at the scale of the metapopulation
            size_t metapopsize = 0u;
            for (size_t p = 0u; p < 2u; ++p) {
                for (auto ind : pops[p].individuals) {

                    const vecDbl traitValues = ind->getTraits();
                    const vecDbl geneticValues = ind->getGeneticValues();
                    for (size_t trait = 0u; trait < 3u; ++trait) {
                        meanPhenotypes[trait][2u] += traitValues[trait];
                        meanGenValues[trait][2u] += geneticValues[trait];
                        pheVariances[trait][2u] += sqr(traitValues[trait]);
                        genVariances[trait][2u] += sqr(geneticValues[trait]);
                    }
                }
                metapopsize += pops[p].getPopSize();
            }
            for (size_t trait = 0u; trait < 2u; ++trait) {
                meanPhenotypes[trait][2u] /= metapopsize;
                pheVariances[trait][2u] /= metapopsize;
                pheVariances[trait][2u] -= sqr(meanPhenotypes[trait][2u]);
                genVariances[trait][2u] /= metapopsize;
                genVariances[trait][2u] -= sqr(meanGenValues[trait][2u]);
            }

            // Assign ecotypes and calculate ecotype-specific means
            for (size_t eco = 0u; eco < 2u; ++eco)
                ecotypes[eco].clear();

            for (size_t p = 0u; p < 2u; ++p) {
                for (auto ind : pops[p].individuals) {
                    size_t group = ind->getEcoTrait() < meanPhenotypes[0u][2u];
                    vecDbl traitValues = ind->getTraits();
                    vecDbl geneticValues = ind->getGeneticValues();
                    for (size_t trait = 0u; trait < 3u; ++trait) {
                        meanPhenotypes[trait][group] += traitValues[trait];
                        meanGenValues[trait][group] += geneticValues[trait];
                        pheVariances[trait][group] += sqr(traitValues[trait]);
                        genVariances[trait][group] += sqr(geneticValues[trait]);
                    }
                    ecotypes[group].push_back(ind);
                }
            }

            for (size_t trait = 0u; trait < 2u; ++trait) {
                for (size_t eco = 0u; eco < 2u; ++eco) {
                    meanPhenotypes[trait][eco] /= ecotypes[eco].size();
                    meanGenValues[trait][eco] /= ecotypes[eco].size();
                    pheVariances[trait][eco] /= ecotypes[eco].size();
                    pheVariances[trait][eco] -= meanPhenotypes[trait][eco];
                    genVariances[trait][eco] /= ecotypes[eco].size();
                    genVariances[trait][eco] -= meanGenValues[trait][eco];
                }
            }

            const size_t n0 = ecotypes[0u].size();
            const size_t n1 = ecotypes[1u].size();
            const vecUns census = { n0, n1, metapopsize };

            // Locus-specific variance decomposition
            for (size_t locus = 0u; locus < genome.nloci; ++locus) {

                double meanAlleleCount = 0.0;
                double varAlleleCount = 0.0;
                vecDbl locusMeanG = zeros(3u);
                vecDbl locusVarG = zeros(3u);
                double covGenValueAlleleCount = 0.0;

                MatUns genotypeCounts = { uzeros(3u), uzeros(3u), uzeros(3u) };
                vecDbl meanGenotypeGenValues = zeros(3u);

                for (size_t eco = 0u; eco < 2u; ++eco) {
                    for (auto ind : ecotypes[eco]) {

                        size_t zyg = ind->getZygosity(locus);
                        double genvalue = ind->getLocusGenValue(locus);

                        meanAlleleCount += zyg;
                        varAlleleCount += sqr(zyg);

                        ++genotypeCounts[2u][zyg];
                        ++genotypeCounts[eco][zyg];
                        meanGenotypeGenValues[zyg] += genvalue;

                        locusMeanG[eco] += genvalue;
                        locusMeanG[2u] += genvalue;
                        locusVarG[eco] += sqr(genvalue);
                        locusVarG[2u] += sqr(genvalue);
                        covGenValueAlleleCount += zyg * genvalue;

                    }
                }

                meanAlleleCount /= metapopsize;
                varAlleleCount /= metapopsize;
                varAlleleCount -= meanAlleleCount;

                for (size_t eco = 0u; eco < 3u; ++eco) {
                    locusMeanG[eco] /= census[eco];
                    locusVarG[eco] /= census[eco];
                    locusVarG[eco] -= sqr(locusMeanG[eco]);
                }

                covGenValueAlleleCount /= metapopsize;
                covGenValueAlleleCount -= meanAlleleCount * locusMeanG[2u];
                for (size_t zyg = 0u; zyg < 3u; ++zyg)
                    meanGenotypeGenValues[zyg] /= genotypeCounts[2u][zyg];

                const size_t trait = genome.traits[locus];

                double avgMutEffect = covGenValueAlleleCount / varAlleleCount;
                vecDbl locusVarA = zeros(3u);
                vecDbl locusMeanBreed = zeros(2u); // 3rd is zero
                locusVarA[2u] = sqr(avgMutEffect) * varAlleleCount;
                addVariances[trait][2u] += locusVarA[2u];

                vecDbl domDeviations = zeros(3u);
                vecDbl breedingValues = zeros(3u);
                vecDbl addExpectations = zeros(3u);
                double locusVarD = 0.0;
                for (size_t zyg = 0u; zyg < 3u; ++zyg) {
                    breedingValues[zyg] = avgMutEffect;
                    breedingValues[zyg] *= (zyg - meanAlleleCount);
                    addExpectations[zyg] = locusMeanG[2u];
                    addExpectations[zyg] -= breedingValues[zyg];
                    domDeviations[zyg] = meanGenotypeGenValues[zyg];
                    domDeviations[zyg] -= addExpectations[zyg];
                    double genotypeSSDeviation = genotypeCounts[2u][zyg];
                    genotypeSSDeviation += sqr(domDeviations[zyg]);
                    locusVarD += genotypeSSDeviation;
                    for (size_t eco = 0u; eco < 2u; ++eco) {
                        size_t n = genotypeCounts[eco][zyg];
                        locusVarA[eco] += n * sqr(breedingValues[zyg]);
                        locusMeanBreed[eco] += n * breedingValues[zyg];
                    }
                }
                for (size_t eco = 0u; eco < 2u; ++eco) {
                    locusMeanBreed[eco] /= census[eco];
                    locusVarA[eco] /= census[eco];
                    locusVarA[eco] -= locusMeanBreed[eco];
                }
                locusVarD /= metapopsize;
                domVariances[trait] += locusVarD;

                double locusVarI = 0.0;
                vecDbl locusMeanDev = zeros(2u);
                vecDbl locusVarN = zeros(3u);
                for (size_t eco = 0u; eco < 2u; ++eco) {
                    for (auto ind : ecotypes[eco]) {

                        const size_t zyg = ind->getZygosity(locus);
                        const double genvalue = ind->getLocusGenValue(locus);

                        double intDeviation = genvalue;
                        intDeviation -= meanGenotypeGenValues[zyg];
                        locusVarI += sqr(intDeviation);

                        double totDeviation = genvalue - addExpectations[zyg];
                        locusVarN[eco] += sqr(totDeviation);
                        locusMeanDev[eco] += totDeviation;

                    }

                    locusVarN[eco] /= ecotypes[eco].size();
                    locusVarN[eco] -= sqr(locusMeanDev[eco]);

                    nadVariances[trait][eco] += locusVarN[eco];

                    double meanBreed = 0.0;
                    double meanBreedSq = 0.0;
                    for (size_t zyg = 0u; zyg < 3u; ++zyg) {
                        size_t n = genotypeCounts[eco][zyg];
                        double brv = breedingValues[zyg];
                        meanBreed += n * brv;
                        meanBreedSq += n * sqr(brv);
                    }
                    meanBreed /= ecotypes[eco].size();
                    meanBreedSq /= ecotypes[eco].size();
                    addVariances[trait][eco] += meanBreedSq - sqr(meanBreed);
                }
                locusVarI /= metapopsize;
                locusVarI *= 2.0;
                intVariances[trait] += locusVarI;

                locusVarN[2u] = locusVarD + 0.5 * locusVarI;
                nadVariances[trait][2u] += locusVarN[2u];

                // Calculate VarS and VarT for genome-wide Fst
                // Go check in a pop gen textbook what these mean
                double Hwithin = 0.0;
                for (size_t eco = 0u; eco < 2u; ++eco) {

                    // Within-ecotype allele frequencies
                    double alleleFreq = genotypeCounts[eco][0u];
                    alleleFreq += 0.5 * genotypeCounts[eco][1u];
                    alleleFreq /= ecotypes[eco].size();
                    varS[trait] += ecotypes[eco].size() * sqr(alleleFreq);
                    double heterozFreq = 2.0 * alleleFreq * (1.0 - alleleFreq);
                    Hwithin += ecotypes[eco].size() * heterozFreq;
                }
                varS[trait] /= metapopsize;
                Hwithin /= metapopsize;

                // Global allele frequency
                double alleleFreq = genotypeCounts[2u][0u];
                alleleFreq += 0.5 * genotypeCounts[2u][1u];
                alleleFreq /= metapopsize;
                varS[trait] -= sqr(alleleFreq);
                varT[trait] += alleleFreq * (1.0 - alleleFreq);

                vecDbl locusVarP = zeros(3u);
                double locusVarE = 1.0;
                for (size_t eco = 0u; eco < 3u; ++eco)
                    locusVarP[eco] = locusVarG[eco] + locusVarE;

                PstScan[locus] = Xst(locusVarP, census);
                GstScan[locus] = Xst(locusVarG, census);
                QstScan[locus] = Xst(locusVarA, census);
                CstScan[locus] = Xst(locusVarN, census);

                // Fst genome scan
                double Htotal = meanAlleleCount * (1.0 - 0.5 * meanAlleleCount);
                FstScan[locus] = 1.0 - Hwithin / Htotal;

            }


            for (size_t trait = 0u; trait < 3u; ++trait) {
                Pst[trait] = Xst(pheVariances[trait], census);
                Gst[trait] = Xst(genVariances[trait], census);
                Qst[trait] = Xst(addVariances[trait], census);
                Cst[trait] = Xst(nadVariances[trait], census);

                Fst[trait] = varS[trait] / varT[trait];
            }

            // Load output to buffer
            loadBuffer(t);

            // Write to files
            for (size_t f = 0u; f < out.names.size(); ++f)
                buffer.write(out.files[f], buffer.fields[f]);


        }

        // Dispersal (if not burnin)
        if (t > 0) {
            Crowd migrants1 = pops[0u].emigrate(dispersal);
            Crowd migrants2 = pops[1u].emigrate(dispersal);
            pops[0u].immigrate(migrants2);
            pops[1u].immigrate(migrants1);
        }

        // Feeding
        if (t > 0) {
            pops[0u].consume();
            pops[1u].consume();
        } else {
            pops[0u].burninConsume();
            pops[1u].burninConsume();
        }


        // Reproduction
        if (t > 0) {
            pops[0u].reproduce(birth, sexsel, genome, networks);
            pops[1u].reproduce(birth, sexsel, genome, networks);
        } else {
            pops[0u].burninReproduce(birth, sexsel, genome, networks, ecosel);
            pops[1u].burninReproduce(birth, sexsel, genome, networks, ecosel);
        }


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
        buffer.add(genVariances[trait][2u]);
        buffer.add(addVariances[trait][2u]);
        buffer.add(domVariances[trait]);
        buffer.add(intVariances[trait]);
        buffer.add(Pst[trait]);
        buffer.add(Gst[trait]);
        buffer.add(Qst[trait]);
        buffer.add(Cst[trait]);
        buffer.add(Fst[trait]);
    }

    buffer.add(getEcoIsolation());
    buffer.add(getSpatialIsolation());
    buffer.add(getMatingIsolation());
}

double MetaPop::getEcoIsolation()
{
    // Ecological isolation is the standard deviation in ecological trait value

    double ei = 0.0;
    const double ecomean = meanPhenotypes[0u][2u];

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










