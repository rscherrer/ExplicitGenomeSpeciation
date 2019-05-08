/*==================================================================================================================================
                                                     analysis.cpp
====================================================================================================================================

C++-code accompanying:	
		 
		(ms. in prep).

Written by:
        G. Sander van Doorn
       	Centre for Ecological and Evolutionary Studies - Theoretical Biology Group
        University of Groningen
        the Netherlands

Program version
		xx/xx/2018	:	

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <queue>
#include "analysis.h"
#include "random.h"
#include "BufferBox.h"
#include "Individual.h"
#include "Population.h"
#include "Genome.h"
#include "Buffer.h"


/*=======================================================================================================
                                     output buffer file
========================================================================================================*/





/*=======================================================================================================
                                     data analysis functions
========================================================================================================*/

double computePostIsolation(const ParameterSet& parameters, const Population& population, const Genome& genome)
{

    // sort out females and males
    std::queue<PInd> females;
    std::vector<PInd> males;
    for(PInd pInd : population.individuals) {
        if(pInd->isFemale(parameters.isFemaleHeteroGamety)) females.push(pInd);
        else males.push_back(pInd);
    }

    const size_t n = males.size();
    if(n == 0u) return 0.0;

    std::array<size_t, 4u> ki {0u, 0u, 0u, 0u};
    std::array<size_t, 4u> kij {0u, 0u, 0u, 0u};

    while(!females.empty())
    {
        PInd fem = females.front();
        females.pop();

        // sample mate at random
        size_t j = rnd::random_int(n);

        size_t f = fem->getEcotype() - 1u;
        size_t m = males[j]->getEcotype() - 1u;

        ++ki[f];
        ++ki[2u + m];

        // create offspring
        PInd offsp = new Individual(fem, males[j], parameters, population, genome);

        // developmental viability
        if(rnd::bernoulli(offsp->getViability())) {
            ++kij[2u * f + m];
        }
    }

    double prod = 1.0;
    for(size_t i = 0u; i < 4u; ++i) {
        if(ki[i] == 0u) return 0.0;
        prod *= ki[i];
    }
    return (1.0 * kij[0u] * kij[3u] - 1.0 * kij[1u] * kij[2u]) / sqrt(prod);

}

double computeMatingIsolation(const ParameterSet& parameters, const Population& population)
{
    // sort out females and males
    std::queue<PInd> females;
    std::vector<PInd> males;
    for(PInd pInd : population.individuals) {
        if(pInd->isFemale(parameters.isFemaleHeteroGamety)) females.push(pInd);
        else males.push_back(pInd);
    }
    
    const size_t n = males.size();
    if(n == 0u) return 0.0;
    
    // allow all females to select a mate from the entire population of males
    std::array<size_t, 4u> ki {0u, 0u, 0u, 0u};
    std::array<size_t, 4u> kij {0u, 0u, 0u, 0u};
    while(!females.empty())
    {
        PInd fem = females.front();
        females.pop();
        
        fem->prepareChoice();
        size_t nMale = rnd::poisson(1.0 / parameters.mateEvaluationCost);
        while(nMale) {
            // assess male
            size_t j  = rnd::random_int(n);
            if(fem->acceptMate(males[j], parameters)) {
                size_t f = fem->getEcotype() - 1u;
                size_t m = males[j]->getEcotype() - 1u;
                ++ki[f];
                ++ki[2u + m];
                ++kij[2u * f + m];
                break;
            }
            --nMale;
        }
    }
    double prod = 1.0;
    for(size_t i = 0u; i < 4u; ++i) {
        if(ki[i] == 0u) return 0.0;
        prod *= ki[i];
    }
    return (1.0 * kij[0u] * kij[3u] - 1.0 * kij[1u] * kij[2u]) / sqrt(prod);
}


void recordData(int t, const std::array<size_t, 7u> &n, const ParameterSet& parameters,
        const std::vector<std::pair<double, double> >& resourceConsumption,
        const std::vector<std::pair<double, double> >& resourceEql,
        std::ofstream& arcFile,
        std::ofstream& datFile,
        std::vector<std::pair<size_t, size_t> >& genderCounts,
        const Population& population,
        const Genome& genome)
{

    // n = (whole pop, hab0, hab1, eco1 hab0, eco2 hab0, eco1 hab1, eco2 hab1)

    // export trait means and sequence to fossil record file
    arcFile << t;
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        arcFile  << '\t' << population.avgG[crctr][0u];
    arcFile << '\t' << to_string(population.individuals.front()->getGenome()) << '\n';
    
    // write output to data file
    size_t nfem = 0u, nmal = 0u;
    for(size_t hab = 0u; hab < parameters.nHabitat; ++hab) {
        nfem += genderCounts[hab].first;
        nmal += genderCounts[hab].second;
    }
    datFile << t
    << '\t' << population.individuals.size()
    << '\t' << nfem
    << '\t' << nmal;
    for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
        datFile << '\t' << genderCounts[hab].first + genderCounts[hab].second;
    for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
        datFile << '\t' << resourceConsumption[hab].first << '\t' << resourceConsumption[hab].second;
    for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
        datFile << '\t' << resourceEql[hab].first << '\t' << resourceEql[hab].second;

    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
        datFile << '\t' << n[1u + 2 * crctr] << '\t' << n[2u + 2 * crctr]
        << '\t' << population.avgG[crctr][0u]
        << '\t' << population.avgG[crctr][1u]
        << '\t' << population.avgG[crctr][2u]
        << '\t' << population.varP[crctr][0u]
        << '\t' << population.varG[crctr][0u]
        << '\t' << population.varA[crctr][0u]
        << '\t' << population.varD[crctr]
        << '\t' << population.varI[crctr][0u]
        << '\t' << population.F_st[crctr]
        << '\t' << population.P_st[crctr]
        << '\t' << population.G_st[crctr]
        << '\t' << population.Q_st[crctr]
        << '\t' << population.C_st[crctr];
    }
    
    double SI, EI, RI, PI;
    const size_t n_1 = n[3u] + n[5u], n_2 = n[4u] + n[6u], n1_ = n[3u] + n[4u], n2_ = n[5u] + n[6u];
    if(n_1 == 0u || n_2 == 0u) SI = EI = RI = PI = 0.0;
    else {
       
        SI = (n1_ == 0u || n2_ == 0u) ? 0.0 : (1.0 * n[3u] * n[6u] - 1.0 * n[4u] * n[5u]) / sqrt(1.0 * n_1 * n_2 * n1_ * n2_);
        EI = population.P_st[0u];
        RI = computeMatingIsolation(parameters, population);
        if(parameters.costIncompat > 0.0) {
            PI = computePostIsolation(parameters, population, genome);
        }
        else {
            PI = 0.0;
        }
    }

    datFile << '\t' << SI << '\t' << EI << '\t' << RI << "\t" << PI;
    datFile << '\n';
    datFile.flush();
    
    // screen output
    std::cout << "t = " << t << ", n = " << population.individuals.size() << ", SI =  " << SI << ", EI = " << EI << ", RI = " << RI << ", PI = " << PI << '\n';
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        std::cout   << "\ttrait." << crctr << " : " << population.avgG[crctr][0u] << " +/- " << sqrt(population.varG[crctr][0u]) << '\n';
}

double Xst(const double &var0, const double &var1, const double &var2, const std::array<size_t, 7u> &n, const double& tiny)
{
    if(var0 < tiny) return 0.0;
    else {
        double varin = (n[1u] * var1 + n[2u] * var2) / n[0u];
        double Xst = 1.0 - varin / var0;
        if(Xst < tiny) Xst = 0.0;
        if(Xst > 1.0 - tiny) Xst = 1.0;
        return Xst;
    }
}

void decomposeVariance(int t,
        const ParameterSet& parameters,
        BufferBox& bufferPointers,
        const Individual::TradeOffPt& breakEvenPoint,
        const std::vector<std::pair<double, double> >& resourceConsumption,
        const std::vector<std::pair<double, double> >& resourceEql,
        std::ofstream& arcFile,
        std::ofstream& datFile,
        std::vector<std::pair<size_t, size_t> >& genderCounts,
        Population& population,
        Genome& genome)
{
    std::array<size_t, 7u> n {population.individuals.size(), 0u, 0u};
    
    // *** genome-wide decomposition of genetic variance ***
    // set initial values
    for(size_t cl = 0u; cl < 3u; ++cl)
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
            population.avgG[crctr][cl] = 0.0;
            population.varP[crctr][cl] = 0.0;
            population.varG[crctr][cl] = 0.0;
        }
    

    // assign ecotype, compute avgG, varG and varP from phenotypic values,
    for(PInd pInd : population.individuals) {
        size_t cl =  pInd->setEcotype(breakEvenPoint);
        ++n[cl];
        ++n[2u + cl + (2u * pInd->getHabitat())];
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
            double g = pInd->getTraitG()[crctr];
            double pp = sqr(pInd->getTraitP()[crctr]);
            population.avgG[crctr][0u] += g;
            population.varG[crctr][0u] += g * g;
            population.varP[crctr][0u] += pp;
            population.avgG[crctr][cl] += g;
            population.varG[crctr][cl] += g * g;
            population.varP[crctr][cl] += pp;
        }
    }
    for(size_t cl = 0u; cl < 3u; ++cl)
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
            if(n[cl] > 1u) {
                double mu = population.avgG[crctr][cl] /= n[cl];
                double aux = (population.varG[crctr][cl] - n[cl] * sqr(mu)) / (n[cl] - 1u);
                population.varG[crctr][cl] = aux > parameters.tiny ? aux : 0.0;
                aux = (population.varP[crctr][cl] / n[cl] - sqr(mu));
                population.varP[crctr][cl] = aux > parameters.tiny ? aux : 0.0;
            }
            else {
                population.avgG[crctr][cl] = population.avgG[crctr][0u];
                population.varP[crctr][cl] = population.varG[crctr][cl] = 0.0;
            }
        }
    
    // *** single locus decomposition of genetic variance ***
    // loop over all loci
    std::vector<double> varE;
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
        size_t nloc = genome.vertices[crctr].size();
        varE[crctr] = parameters.nLoci ? sqr(parameters.scaleE[crctr]) / nloc : 0.0;
    }
    
    for(size_t i = 0u; i < parameters.nLoci; ++i) {

        // determine genotype frequencies and genetic values
        std::array<size_t, 3u> m {0u, 0u, 0u};
        std::array<double, 3u> dev {0.0, 0.0, 0.0};
        std::array<double, 3u> sumu;
        for(size_t cl = 0u; cl < 3u; ++cl) {
            genome.characterLocus[i].meanEffect[cl] = 0.0;
            genome.characterLocus[i].varG[cl] = 0.0;
            sumu[cl] = 0u;
        }
        for(PInd pInd : population.individuals) {
            size_t u = pInd->getTraitLocus()[i].alleleCount;
            double g = pInd->getTraitLocus()[i].geneticValue;
            size_t cl = pInd->getEcotype();
            ++m[u];
            dev[u] += g;
            sumu[0u] += u;
            sumu[cl] += u;
            genome.characterLocus[i].meanEffect[0u] += g;
            genome.characterLocus[i].meanEffect[cl] += g;
            genome.characterLocus[i].varG[0u] += g * g;
            genome.characterLocus[i].varG[cl] += g * g;
        }
        // allele frequency, mean effect and genetic variance
        size_t ci = genome.characterLocus[i].character;
        double varEi = varE[ci];
        for(size_t cl = 0u; cl < 3u; ++cl) {
            if(n[cl] > 1u) {
                double mu = genome.characterLocus[i].meanEffect[cl] /= n[cl];
                double aux = (genome.characterLocus[i].varG[cl] - n[cl] * sqr(mu)) / (n[cl] - 1u);
                genome.characterLocus[i].varG[cl] = aux > parameters.tiny ? aux : 0.0;
                sumu[cl] /= n[cl];
                aux = 0.5 * sumu[cl];
                if(aux < parameters.tiny) aux = 0.0;
                if(aux > 1.0 - parameters.tiny) aux = 1.0;
                genome.characterLocus[i].alleleFrequency[cl] = aux;
            }
            else {
                genome.characterLocus[i].meanEffect[cl] = genome.characterLocus[i].meanEffect[0u];
                genome.characterLocus[i].varG[cl] = 0.0;
                genome.characterLocus[i].alleleFrequency[cl] = genome.characterLocus[i].alleleFrequency[0u];
            }
            genome.characterLocus[i].varP[cl] = genome.characterLocus[i].varG[cl] + varEi;
        }
        double avgu = sumu[0u];

        // linear regression of genetic value on allele count
        double sumuu = 0.0, sumug = 0.0;
        std::array<double, 3u> frq;
        for(size_t u = 0u; u < 3u; ++u){
            sumuu += m[u] * u * u;
            sumug += dev[u] * u;
            if(m[u]) dev[u] /= m[u];
            frq[u] = m[u] * 1.0 / n[0u];
        }

        // variance in allele count (= 2pq in Hardy Weinberg Eql)
        double mu = genome.characterLocus[i].meanEffect[0u];
        double varuu = (sumuu - n[0u] * sqr(avgu)) / (n[0] - 1u);

        // covariance between genetic value and allele count
        double covug = (sumug - n[0u] * mu * avgu) / (n[0] - 1u);

        // average effect
        double alpha = genome.characterLocus[i].avgEffectOfSubstitution =
            varuu > parameters.tiny ? covug / varuu : 0.0;

        // additive genetic variance
        genome.characterLocus[i].varA[0u] = sqr(alpha) * varuu;

        // observed and expected heterozygosity ~ F_it
        double H_t = avgu * (1.0 - 0.5 * avgu);
        double Fit = H_t > parameters.tiny ? 1.0 - frq[1u] / H_t : 0.0;
        if(Fit < parameters.tiny) Fit = 0.0;
        if(Fit > 1.0 - parameters.tiny) Fit = 1.0;
        genome.characterLocus[i].F_it = Fit;

        // mean heterozygosity in subpopulations
        double H_within = 0.0;
        for(size_t cl = 1u; cl < 3u; ++cl) {
            double pi = genome.characterLocus[i].alleleFrequency[cl];
            H_within += n[cl] * 2.0 * pi * (1.0 - pi);
        }
        H_within /= n[0u];

        // F_st and F_is
        double Fis, Fst = H_t > parameters.tiny ? 1.0 - H_within / H_t : 0.0;
        if(Fst < parameters.tiny) Fst = 0.0;
        if(Fst > 1.0 - parameters.tiny) {
            Fst = 1.0;
            Fis = 0.0;
        }
        else Fis = 1.0 - (1.0 - Fit) / (1.0 - Fst);
        if(Fis < parameters.tiny) Fis = 0.0;
        if(Fis > 1.0 - parameters.tiny) Fis = 1.0;
        genome.characterLocus[i].F_is = Fis;
        genome.characterLocus[i].F_st = Fst;
        
        // dominance deviations and dominance variance
        genome.characterLocus[i].varD = 0.0;
        for(size_t u = 0u; u < 3u; ++u) {
            dev[u] -= mu + alpha * (u - avgu);
            genome.characterLocus[i].varD += frq[u] * sqr(dev[u]);
        }

        // epistatic deviations and epistatic variance
        // additive genetic variation within clusters
        double sumeps2 = 0.0;
        std::array<double, 3u> sumbrv1 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumbrv2 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumdva1 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumdva2 {0.0, 0.0, 0.0};
        for(PInd pInd : population.individuals) {
            size_t cl = pInd->getEcotype();
            size_t u = pInd->getTraitLocus()[i].alleleCount;
            double g = pInd->getTraitLocus()[i].geneticValue;
            double brv = alpha * (u - avgu);    // breeding value
            double dva = g - (mu + brv);        // deviation from additivity
            double eps = dva - dev[u];          // deviation by epistasis
            sumbrv1[cl] += brv;
            sumbrv2[cl] += brv * brv;
            sumdva1[cl] += dva;
            sumdva2[cl] += dva * dva;
            sumeps2 += eps * eps;
        }
        genome.characterLocus[i].varI[0u] = 2.0 * sumeps2 / (n[0] - 1u);
        for(size_t cl = 1u; cl < 3u; ++cl) {
            if(n[cl] > 1u) {
                sumbrv1[cl] /= n[cl];
                sumdva1[cl] /= n[cl];
                double aux = (sumbrv2[cl] - n[cl] * sqr(sumbrv1[cl])) / (n[cl] - 1u);
                genome.characterLocus[i].varA[cl] = aux > parameters.tiny ? aux : 0.0;
                aux = (sumdva2[cl] - n[cl] * sqr(sumdva1[cl])) / (n[cl] - 1u);
                genome.characterLocus[i].varI[cl] = aux > parameters.tiny ? aux : 0.0;
            }
            else {
                genome.characterLocus[i].varA[cl] = 0.0;
                genome.characterLocus[i].varI[cl] = 0.0;
            }
        }
        
        // Pst, Gst, Qst and Cst
        genome.characterLocus[i].P_st =
        Xst(genome.characterLocus[i].varP[0u],
            genome.characterLocus[i].varP[1u],
            genome.characterLocus[i].varP[2u], n, parameters.tiny);
        genome.characterLocus[i].G_st =
            Xst(genome.characterLocus[i].varG[0u],
                genome.characterLocus[i].varG[1u],
                genome.characterLocus[i].varG[2u], n, parameters.tiny);
        genome.characterLocus[i].Q_st =
            Xst(genome.characterLocus[i].varA[0u],
                genome.characterLocus[i].varA[1u],
                genome.characterLocus[i].varA[2u], n, parameters.tiny);
        genome.characterLocus[i].C_st =
            Xst(genome.characterLocus[i].varD + 0.5 * genome.characterLocus[i].varI[0u],
                genome.characterLocus[i].varI[1u],
                genome.characterLocus[i].varI[2u], n, parameters.tiny);
    }

    // write data to output buffer
    for(size_t i = 0u; i < parameters.nLoci; ++i) {
        (*bufferPointers.bufferFreq)[i] = genome.characterLocus[i].alleleFrequency[0u];
        (*bufferPointers.bufferF_it)[i] = genome.characterLocus[i].F_it;
        (*bufferPointers.bufferF_is)[i] = genome.characterLocus[i].F_is;
        (*bufferPointers.bufferF_st)[i] = genome.characterLocus[i].F_st;
        (*bufferPointers.bufferP_st)[i] = genome.characterLocus[i].P_st;
        (*bufferPointers.bufferG_st)[i] = genome.characterLocus[i].G_st;
        (*bufferPointers.bufferQ_st)[i] = genome.characterLocus[i].Q_st;
        (*bufferPointers.bufferC_st)[i] = genome.characterLocus[i].C_st;
        (*bufferPointers.bufferVarP)[i] = genome.characterLocus[i].varP[0u];
        (*bufferPointers.bufferVarG)[i] = genome.characterLocus[i].varG[0u];
        (*bufferPointers.bufferVarA)[i] = genome.characterLocus[i].varA[0u];
        (*bufferPointers.bufferVarD)[i] = genome.characterLocus[i].varD;
        (*bufferPointers.bufferVarI)[i] = genome.characterLocus[i].varI[0u];
    }
    bufferPointers.bufferFreq->flush(parameters);
    bufferPointers.bufferF_it->flush(parameters);
    bufferPointers.bufferF_is->flush(parameters);
    bufferPointers.bufferF_st->flush(parameters);
    bufferPointers.bufferP_st->flush(parameters);
    bufferPointers.bufferG_st->flush(parameters);
    bufferPointers.bufferQ_st->flush(parameters);
    bufferPointers.bufferC_st->flush(parameters);
    bufferPointers.bufferVarP->flush(parameters);
    bufferPointers.bufferVarG->flush(parameters);
    bufferPointers.bufferVarA->flush(parameters);
    bufferPointers.bufferVarD->flush(parameters);
    bufferPointers.bufferVarI->flush(parameters);

    // *** genome-wide decomposition of genetic variance (continued) ***
    // compute varA, varD, varI, and Fst by accumulating single locus contributions
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
        for(size_t cl = 0u; cl < 3u; ++cl) {
            population.varA[crctr][cl] = 0.0;
            population.varI[crctr][cl] = 0.0;
        }
        population.varD[crctr] = 0.0;
        double var_t = 0.0, var_s = 0.0;
        for(size_t i : genome.vertices[crctr]) {
            for(size_t cl = 0u; cl < 3u; ++cl) {
                population.varA[crctr][cl] += genome.characterLocus[i].varA[cl];
                population.varI[crctr][cl] += genome.characterLocus[i].varI[cl];
            }
            population.varD[crctr] += genome.characterLocus[i].varD;
            double p0 = genome.characterLocus[i].alleleFrequency[0u];
            double p1 = genome.characterLocus[i].alleleFrequency[1u];
            double p2 = genome.characterLocus[i].alleleFrequency[2u];
            var_s += (n[1u] * sqr(p1) + n[2u] * sqr(p2)) / n[0u] - sqr(p0);
            var_t += p0 * (1.0 - p0);
        }
        population.F_st[crctr] = var_t < parameters.tiny ? 0.0 : var_s / var_t;
        population.P_st[crctr] = Xst(population.varP[crctr][0u],
                                      population.varP[crctr][1u],
                                      population.varP[crctr][2u], n, parameters.tiny);
        population.G_st[crctr] = Xst(population.varG[crctr][0u],
                                      population.varG[crctr][1u],
                                      population.varG[crctr][2u], n, parameters.tiny);
        population.Q_st[crctr] = Xst(population.varA[crctr][0u],
                                      population.varA[crctr][1u],
                                      population.varA[crctr][2u], n, parameters.tiny);
        population.C_st[crctr] =
            Xst(population.varD[crctr] + 0.5 * population.varI[crctr][0u],
                population.varI[crctr][1u],
                population.varI[crctr][2u], n, parameters.tiny);
    }
    recordData(t, n, parameters, resourceConsumption, resourceEql, arcFile, datFile, genderCounts, population, genome);
}

void analyseNetwork(int t, const ParameterSet& parameters, const Population& population, const Genome& genome)
{
    const char sep = ',';
    // *** node properties ***
    // open network link file
    std::ostringstream oss1;
    oss1 << "simulation_" << parameters.seed << '_' << t << "_nodes.csv";
    std::ofstream ofs(oss1.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    
    ofs << "id" << sep << "character" << sep << "linkage.group" << sep
        << "degree" << sep << "location" << sep
        << "effect.size" << sep << "dominance.coefficient" << sep
        << "allele.freq" << sep << "mean.effect" << sep
        << "avg.effect.substitution" << sep << "varP" << sep
        << "varG" << sep << "varA" << sep << "varD" << sep
        << "varI" << sep << "exp.heterozygosity" << sep
        << "Fit" << sep << "Fis" << sep << "Fst" << sep
        << "Gst" << sep << "Qst" << sep << "Cst" << '\n';
    
    for(size_t i = 0u; i < parameters.nLoci; ++i) {
        const size_t crctr = genome.characterLocus[i].character;
        double pi = genome.characterLocus[i].alleleFrequency[0u];
        ofs << i << sep
            << genome.characterLocus[i].character << sep
            << genome.characterLocus[i].linkageGroup << sep
            << genome.characterLocus[i].edges.size() << sep
            << genome.characterLocus[i].location << sep
            << parameters.scaleA[crctr] * genome.characterLocus[i].effectSize << sep
            << parameters.scaleD[crctr] * genome.characterLocus[i].dominanceCoeff << sep
            << pi << sep
            << genome.characterLocus[i].meanEffect[0u] << sep
            << genome.characterLocus[i].avgEffectOfSubstitution << sep
            << genome.characterLocus[i].varP[0u] << sep
            << genome.characterLocus[i].varG[0u] << sep
            << genome.characterLocus[i].varA[0u] << sep
            << genome.characterLocus[i].varD << sep
            << genome.characterLocus[i].varI[0u] << sep
            << 2.0 * pi * (1.0 - pi) << sep
            << genome.characterLocus[i].F_it << sep
            << genome.characterLocus[i].F_is << sep
            << genome.characterLocus[i].F_st << sep
            << genome.characterLocus[i].G_st << sep
            << genome.characterLocus[i].Q_st << sep
            << genome.characterLocus[i].C_st <<'\n';
    }
    ofs.close();

    // *** edge properties ***
    // open network link file
    std::ostringstream oss2;
    oss2 << "simulation_" << parameters.seed << '_' << t << "_links.csv";
    ofs.open(oss2.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    ofs << "from" << sep << "to" << sep << "effect.size" << sep
        << "interaction.weight" << sep
        << "gene.freq.correlation" << sep
        << "genetic.correlation" << sep
        << "additive.genetic.correlation" << '\n';
    
    // determine edge properties
    const size_t n = population.individuals.size();
    for(size_t i = 0u; i < parameters.nLoci; ++i) {
        const size_t crctr = genome.characterLocus[i].character;
        for(const std::pair<size_t, double> &edge : genome.characterLocus[i].edges) {
            const size_t j = edge.first;
            const double eij = parameters.scaleI[crctr] * edge.second;
            double sumpipj = 0.0, sumgigj = 0.0, sumbibj = 0.0,
                sumxi = 0.0, sumxj = 0.0, sumxixi = 0.0, sumxjxj;
            double pi = genome.characterLocus[i].alleleFrequency[0u];
            double pj = genome.characterLocus[j].alleleFrequency[0u];
            double alphai = genome.characterLocus[i].avgEffectOfSubstitution;
            double alphaj = genome.characterLocus[j].avgEffectOfSubstitution;
            for(PInd pInd : population.individuals) {
                size_t ui = pInd->getTraitLocus()[i].alleleCount;
                size_t uj = pInd->getTraitLocus()[j].alleleCount;
                double xi = pInd->getTraitLocus()[i].expression;
                double xj = pInd->getTraitLocus()[j].expression;
                double gi = pInd->getTraitLocus()[i].geneticValue;
                double gj = pInd->getTraitLocus()[j].geneticValue;
                double bi = (ui - 2.0 * pi) * alphai;
                double bj = (uj - 2.0 * pj) * alphaj;
                sumpipj += 0.25 * ui * uj;
                sumgigj += gi * gj;
                sumbibj += bi * bj;
                sumxi += xi;
                sumxj += xj;
                sumxixi += xi * xi;
                sumxjxj += xj * xj;
            }
            
            // variation in average effect due to epistasis
            sumxi /= n;
            sumxj /= n;
            sumxixi = (sumxixi - n * sumxi * sumxi) / (n - 1u);
            sumxjxj = (sumxjxj - n * sumxj * sumxj) / (n - 1u);
            double ai = parameters.scaleA[crctr] * genome.characterLocus[i].effectSize;
            double aj = parameters.scaleA[crctr] * genome.characterLocus[j].effectSize;
            double weightij = fabs(aj) < parameters.tiny ? 0.0 :
                sqrt(sqr(genome.characterLocus[j].avgEffectOfSubstitution * eij / aj)
                    * sumxixi);
            double weightji = fabs(ai) < parameters.tiny ? 0.0 :
                sqrt(sqr(genome.characterLocus[i].avgEffectOfSubstitution * eij / ai)
                    * sumxjxj);
        
            // correlation between allele frequencies ~ population LD
            sumpipj = (sumpipj - n * pi * pj) / (n - 1u);
            double norm = sqrt(pi * (1.0 - pi) * pj * (1.0 - pj));
            double rp = norm > parameters.tiny ? sumpipj / norm : 0.0;
            
            // correlation between genetic values ~ genetic covariance
            double mui = genome.characterLocus[i].meanEffect[0u];
            double muj = genome.characterLocus[j].meanEffect[0u];
            double vari = genome.characterLocus[i].varG[0u];
            double varj = genome.characterLocus[j].varG[0u];
            sumgigj = (sumgigj - n * mui * muj) / (n - 1u);
            norm = sqrt(vari * varj);
            double rG = norm > parameters.tiny ? sumgigj / norm : 0.0;
            
            // correlation between breeding values ~ additive genetic covariance
            vari = genome.characterLocus[i].varA[0u];
            varj = genome.characterLocus[j].varA[0u];
            sumbibj /= (n - 1u);
            norm = sqrt(vari * varj);
            double rA = norm > parameters.tiny ? sumbibj / norm : 0.0;
            
            // write edge data to file
            ofs << i << sep << j << sep << eij << sep << weightij << sep
                << rp << sep << rG << sep << rA << '\n';
            ofs << j << sep << i << sep << eij << sep << weightji << sep
                << rp << sep << rG << sep << rA << '\n';
        }
    }
    ofs.close();
    
    // *** trait distributions ***
    // open file
    std::ostringstream oss3;
    oss3 << "simulation_" << parameters.seed << '_' << t << "_dump.csv";
    ofs.open(oss3.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    ofs << sep << "ecotype" << sep << "habitat" << sep << "attack.rate.1" << sep << "attack.rate.2";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << "trait." << crctr;
    ofs << '\n';
    size_t i = 0u;
    for(PInd pInd : population.individuals) {
        ofs << i  << sep << pInd->getEcotype() << sep << pInd->getHabitat()
            << sep << pInd->getAttackRate().first << sep << pInd->getAttackRate().second;
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
            ofs << sep << pInd->getTraitP()[crctr];
        ofs << '\n';
        ++i;
    }
    ofs.close();
}