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

extern int tSavDat, tGetDat, tBurnIn;
extern unsigned int seed;
extern double mateEvaluationCost;
extern std::array<double, nCharacter> scaleA, scaleD, scaleI, scaleE;
extern std::list<PInd> population;
extern std::ofstream datFile, arcFile;
extern Buffer *bufferFreq, *bufferF_it, *bufferF_is, *bufferF_st,
    *bufferP_st, *bufferG_st, *bufferQ_st, *bufferC_st,
    *bufferVarP, *bufferVarG, *bufferVarA, *bufferVarD, *bufferVarI, *bufferCovL;
extern std::array<std::pair<double, double>, nHabitat> resourceConsumption, resourceEql;
extern std::array<std::pair<size_t, size_t>, nHabitat> genderCounts;
extern Individual::TradeOffPt breakEvenPoint;

/*=======================================================================================================
                                     output buffer file
========================================================================================================*/

const char Buffer::sep = ',';

Buffer::Buffer(const std::string &str) :
    i(0u), k(0u), t(-tBurnIn), n(static_cast<size_t>(tSavDat / tGetDat)), label(str)
{
    data = std::vector< std::vector<double> >(n, std::vector<double>(nLoci, 0.0));
    
    // open table data file
    std::ostringstream oss;
    oss << "simulation_" << seed << "_table_" << label << ".csv";
    ofs.open(oss.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in Buffer::Buffer()");
    
    // write header
    ofs << sep << "generation";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << "loc." << j;
    ofs << '\n';
    ofs << "character" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].character;
    ofs << '\n';
    ofs << "linkage.group" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].linkageGroup;
    ofs << '\n';
    ofs << "degree" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].edges.size();
    ofs << '\n';
    ofs << "location" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].location;
    ofs << '\n';
    ofs << "effect.size" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].effectSize;
    ofs << '\n';
    ofs << "dominance.coeff" << sep << "NA";
    for(size_t j = 0u; j < nLoci; ++j)
        ofs << sep << Individual::characterLocus[j].dominanceCoeff;
    ofs << '\n';
}

void Buffer::flush()
{
    t += tGetDat;
    ++i;
    if(t % tSavDat == 0u) {
        ++k;
        ofs << "point." << k << sep << t - (tSavDat >> 1u);
        for(size_t j = 0u; j < nLoci; ++j) {
            double sum = 0.0;
            for(i = 0u; i < n; ++i) sum += data[i][j];
            ofs << sep << sum / n;
        }
        ofs << '\n';
        i = 0u;
    }
    ofs.flush();
}

/*=======================================================================================================
                                     data analysis functions
========================================================================================================*/

double computeMatingIsolation()
{
    // sort out females and males
    std::queue<PInd> females;
    std::vector<PInd> males;
    for(PInd pInd : population) {
        if(pInd->isFemale()) females.push(pInd);
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
        size_t nMale = rnd::poisson(1.0 / mateEvaluationCost);
        while(nMale) {
            // assess male
            size_t j  = rnd::random_int(n);
            if(fem->acceptMate(males[j])) {
                size_t f = fem->ecotype - 1u;
                size_t m = males[j]->ecotype - 1u;
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


void recordData(int t, const std::array<size_t, 7u> &n)
{
    // export trait means and sequence to fossil record file
    arcFile << t;
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        arcFile  << '\t' << Individual::avgG[crctr][0u];
    arcFile << '\t' << population.front()->genome.to_string() << '\n';
    
    // write output to data file
    size_t nfem = 0u, nmal = 0u;
    for(size_t hab = 0u; hab < nHabitat; ++hab) {
        nfem += genderCounts[hab].first;
        nmal += genderCounts[hab].second;
    }
    datFile << t
    << '\t' << population.size()
    << '\t' << nfem
    << '\t' << nmal;
    for(size_t hab = 0u; hab < nHabitat; ++hab)
        datFile << '\t' << genderCounts[hab].first + genderCounts[hab].second;
    for(size_t hab = 0u; hab < nHabitat; ++hab)
        datFile << '\t' << resourceConsumption[hab].first << '\t' << resourceConsumption[hab].second;
    for(size_t hab = 0u; hab < nHabitat; ++hab)
        datFile << '\t' << resourceEql[hab].first << '\t' << resourceEql[hab].second;
    
    
    datFile << '\t' << n[1u] << '\t' << n[2u];
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
        datFile << '\t' << Individual::avgG[crctr][0u]
        << '\t' << Individual::avgG[crctr][1u]
        << '\t' << Individual::avgG[crctr][2u]
        << '\t' << Individual::varP[crctr][0u]
        << '\t' << Individual::varG[crctr][0u]
        << '\t' << Individual::varA[crctr][0u]
        << '\t' << Individual::varD[crctr]
        << '\t' << Individual::varI[crctr][0u]
        << '\t' << Individual::covL[crctr]
        << '\t' << Individual::F_st[crctr]
        << '\t' << Individual::P_st[crctr]
        << '\t' << Individual::G_st[crctr]
        << '\t' << Individual::Q_st[crctr]
        << '\t' << Individual::C_st[crctr];
    }
    
    double SI, EI, RI;
    const size_t n_1 = n[3u] + n[5u], n_2 = n[4u] + n[6u], n1_ = n[3u] + n[4u], n2_ = n[5u] + n[6u];
    if(n_1 == 0u || n_2 == 0u) SI = EI = RI = 0.0;
    else {
       
        SI = (n1_ == 0u || n2_ == 0u) ? 0.0 : (1.0 * n[3u] * n[6u] - 1.0 * n[4u] * n[5u]) / sqrt(1.0 * n_1 * n_2 * n1_ * n2_);
        EI = Individual::P_st[0u];
        RI = computeMatingIsolation();
    }
    datFile << '\t' << SI << '\t' << EI << '\t' << RI;
    datFile << '\n';
    datFile.flush();
    
    // screen output
    std::cout << "t = " << t << ", n = " << population.size() << ", SI =  " << SI << ", EI = " << EI << ", RI = " << RI << '\n';
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        std::cout   << "\ttrait." << crctr << " : " << Individual::avgG[crctr][0u] << " +/- " << sqrt(Individual::varG[crctr][0u]) << '\n';
}

double Xst(const double &var0, const double &var1, const double &var2, const std::array<size_t, 7u> &n)
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

void decomposeVariance(int t)
{
    std::array<size_t, 7u> n {population.size(), 0u, 0u};
    
    // *** genome-wide decomposition of genetic variance ***
    // set initial values
    for(size_t cl = 0u; cl < 3u; ++cl)
        for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
            Individual::avgG[crctr][cl] = 0.0;
            Individual::varP[crctr][cl] = 0.0;
            Individual::varG[crctr][cl] = 0.0;
        }
    

    // assign ecotype, compute avgG, varG and varP from phenotypic values,
    for(PInd pInd : population) {
        size_t cl =  pInd->setEcotype(breakEvenPoint);
        ++n[cl];
        ++n[2u + cl + (2u * pInd->habitat)];
        for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
            double g = pInd->traitG[crctr];
            double pp = sqr(pInd->traitP[crctr]);
            Individual::avgG[crctr][0u] += g;
            Individual::varG[crctr][0u] += g * g;
            Individual::varP[crctr][0u] += pp;
            Individual::avgG[crctr][cl] += g;
            Individual::varG[crctr][cl] += g * g;
            Individual::varP[crctr][cl] += pp;
        }
    }
    for(size_t cl = 0u; cl < 3u; ++cl)
        for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
            if(n[cl] > 1u) {
                double mu = Individual::avgG[crctr][cl] /= n[cl];
                double aux = (Individual::varG[crctr][cl] - n[cl] * sqr(mu)) / (n[cl] - 1u);
                Individual::varG[crctr][cl] = aux > tiny ? aux : 0.0;
                aux = (Individual::varP[crctr][cl] / n[cl] - sqr(mu));
                Individual::varP[crctr][cl] = aux > tiny ? aux : 0.0;
            }
            else {
                Individual::avgG[crctr][cl] = Individual::avgG[crctr][0u];
                Individual::varP[crctr][cl] = Individual::varG[crctr][cl] = 0.0;
            }
        }
    
    // *** single locus decomposition of genetic variance ***
    // loop over all loci
    std::array<double, nCharacter> varE;
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
        size_t nloc = Individual::vertices[crctr].size();
        varE[crctr] = nLoci ? sqr(scaleE[crctr]) / nloc : 0.0;
    }
    
    for(size_t i = 0u; i < nLoci; ++i) {
        // determine genotype frequencies and genetic values
        std::array<size_t, 3u> m {0u, 0u, 0u};
        std::array<double, 3u> dev {0.0, 0.0, 0.0};
        std::array<double, 3u> sumu;
        for(size_t cl = 0u; cl < 3u; ++cl) {
            Individual::characterLocus[i].meanEffect[cl] = 0.0;
            Individual::characterLocus[i].varG[cl] = 0.0;
            sumu[cl] = 0u;
        }
        for(PInd pInd : population) {
            size_t u = pInd->traitLocus[i].alleleCount;
            double g = pInd->traitLocus[i].geneticValue;
            size_t cl = pInd->ecotype;
            ++m[u];
            dev[u] += g;
            sumu[0u] += u;
            sumu[cl] += u;
            Individual::characterLocus[i].meanEffect[0u] += g;
            Individual::characterLocus[i].meanEffect[cl] += g;
            Individual::characterLocus[i].varG[0u] += g * g;
            Individual::characterLocus[i].varG[cl] += g * g;
        }
        // allele frequency, mean effect and genetic variance
        size_t ci = Individual::characterLocus[i].character;
        double varEi = varE[ci];
        for(size_t cl = 0u; cl < 3u; ++cl) {
            if(n[cl] > 1u) {
                double mu = Individual::characterLocus[i].meanEffect[cl] /= n[cl];
                double aux = (Individual::characterLocus[i].varG[cl] - n[cl] * sqr(mu)) / (n[cl] - 1u);
                Individual::characterLocus[i].varG[cl] = aux > tiny ? aux : 0.0;
                sumu[cl] /= n[cl];
                aux = 0.5 * sumu[cl];
                if(aux < tiny) aux = 0.0;
                if(aux > 1.0 - tiny) aux = 1.0;
                Individual::characterLocus[i].alleleFrequency[cl] = aux;
            }
            else {
                Individual::characterLocus[i].meanEffect[cl] = Individual::characterLocus[i].meanEffect[0u];
                Individual::characterLocus[i].varG[cl] = 0.0;
                Individual::characterLocus[i].alleleFrequency[cl] = Individual::characterLocus[i].alleleFrequency[0u];
            }
            Individual::characterLocus[i].varP[cl] = Individual::characterLocus[i].varG[cl] + varEi;
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
        double mu = Individual::characterLocus[i].meanEffect[0u];
        double varuu = (sumuu - n[0u] * sqr(avgu)) / (n[0] - 1u);
        // covariance between genetic value and allele count
        double covug = (sumug - n[0u] * mu * avgu) / (n[0] - 1u);
        // average effect
        double alpha = Individual::characterLocus[i].avgEffectOfSubstitution =
            varuu > tiny ? covug / varuu : 0.0;
        // additive genetic variance
        Individual::characterLocus[i].varA[0u] = sqr(alpha) * varuu;
        // observed and expected heterozygosity ~ F_it
        double H_t = avgu * (1.0 - 0.5 * avgu);
        double Fit = H_t > tiny ? 1.0 - frq[1u] / H_t : 0.0;
        if(Fit < tiny) Fit = 0.0;
        if(Fit > 1.0 - tiny) Fit = 1.0;
        Individual::characterLocus[i].F_it = Fit;
        // mean heterozygosity in subpopulations
        double H_within = 0.0;
        for(size_t cl = 1u; cl < 3u; ++cl) {
            double pi = Individual::characterLocus[i].alleleFrequency[cl];
            H_within += n[cl] * 2.0 * pi * (1.0 - pi);
        }
        H_within /= n[0u];
        // F_st and F_is
        double Fis, Fst = H_t > tiny ? 1.0 - H_within / H_t : 0.0;
        if(Fst < tiny) Fst = 0.0;
        if(Fst > 1.0 - tiny) {
            Fst = 1.0;
            Fis = 0.0;
        }
        else Fis = 1.0 - (1.0 - Fit) / (1.0 - Fst);
        if(Fis < tiny) Fis = 0.0;
        if(Fis > 1.0 - tiny) Fis = 1.0;
        Individual::characterLocus[i].F_is = Fis;
        Individual::characterLocus[i].F_st = Fst;
        
        // dominance deviations and dominance variance
        Individual::characterLocus[i].varD = 0.0;
        for(size_t u = 0u; u < 3u; ++u) {
            dev[u] -= mu + alpha * (u - avgu);
            Individual::characterLocus[i].varD += frq[u] * sqr(dev[u]);
        }
        // epistatic deviations and epistatic variance
        // additive genetic variation within clusters
        double sumeps2 = 0.0;
        std::array<double, 3u> sumbrv1 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumbrv2 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumdva1 {0.0, 0.0, 0.0};
        std::array<double, 3u> sumdva2 {0.0, 0.0, 0.0};
        for(PInd pInd : population) {
            size_t cl = pInd->ecotype;
            size_t u = pInd->traitLocus[i].alleleCount;
            double g = pInd->traitLocus[i].geneticValue;
            double brv = alpha * (u - avgu);    // breeding value
            double dva = g - (mu + brv);        // deviation from additivity
            double eps = dva - dev[u];          // deviation by epitasis
            sumbrv1[cl] += brv;
            sumbrv2[cl] += brv * brv;
            sumdva1[cl] += dva;
            sumdva2[cl] += dva * dva;
            sumeps2 += eps * eps;
        }
        Individual::characterLocus[i].varI[0u] = 2.0 * sumeps2 / (n[0] - 1u);
        for(size_t cl = 1u; cl < 3u; ++cl) {
            if(n[cl] > 1u) {
                sumbrv1[cl] /= n[cl];
                sumdva1[cl] /= n[cl];
                double aux = (sumbrv2[cl] - n[cl] * sqr(sumbrv1[cl])) / (n[cl] - 1u);
                Individual::characterLocus[i].varA[cl] = aux > tiny ? aux : 0.0;
                aux = (sumdva2[cl] - n[cl] * sqr(sumdva1[cl])) / (n[cl] - 1u);
                Individual::characterLocus[i].varI[cl] = aux > tiny ? aux : 0.0;
            }
            else {
                Individual::characterLocus[i].varA[cl] = 0.0;
                Individual::characterLocus[i].varI[cl] = 0.0;
            }
        }
        
        // Pst, Gst, Qst and Cst
        Individual::characterLocus[i].P_st =
        Xst(Individual::characterLocus[i].varP[0u],
            Individual::characterLocus[i].varP[1u],
            Individual::characterLocus[i].varP[2u], n);
        Individual::characterLocus[i].G_st =
            Xst(Individual::characterLocus[i].varG[0u],
                Individual::characterLocus[i].varG[1u],
                Individual::characterLocus[i].varG[2u], n);
        Individual::characterLocus[i].Q_st =
            Xst(Individual::characterLocus[i].varA[0u],
                Individual::characterLocus[i].varA[1u],
                Individual::characterLocus[i].varA[2u], n);
        Individual::characterLocus[i].C_st =
            Xst(Individual::characterLocus[i].varD + 0.5 * Individual::characterLocus[i].varI[0u],
            Individual::characterLocus[i].varI[1u],
            Individual::characterLocus[i].varI[2u], n);
    }
    
    // linkage deviations and genetic variance due to linkage
    for(size_t i = 0u; i < nLoci; ++i) {
        size_t crctr = Individual::characterLocus[i].character;
        double covGij = 0.0;
        for(PInd pInd : population) {
            double gi = pInd->traitLocus[i].geneticValue;
            double sumGj = 0.0;
            for(size_t j : Individual::vertices[crctr]) {
                if(i == j) continue;
                sumGj += pInd->traitLocus[j].geneticValue;
            }
            covGij += gi * sumGj;
        }
        double mui = Individual::characterLocus[i].meanEffect[0u], mujSum = 0.0;
        for(size_t j : Individual::vertices[crctr]) {
            if(i == j) continue;
            mujSum += Individual::characterLocus[j].meanEffect[0u];
        }
        covGij = (covGij - n[0u] * mui * mujSum) / (n[0u] - 1u);
        Individual::characterLocus[i].covL = covGij - 0.5 * Individual::characterLocus[i].varI[0u]; // can be negative
    }
    // write data to output buffer
    for(size_t i = 0u; i < nLoci; ++i) {
        (*bufferFreq)[i] = Individual::characterLocus[i].alleleFrequency[0u];
        (*bufferF_it)[i] = Individual::characterLocus[i].F_it;
        (*bufferF_is)[i] = Individual::characterLocus[i].F_is;
        (*bufferF_st)[i] = Individual::characterLocus[i].F_st;
        (*bufferP_st)[i] = Individual::characterLocus[i].P_st;
        (*bufferG_st)[i] = Individual::characterLocus[i].G_st;
        (*bufferQ_st)[i] = Individual::characterLocus[i].Q_st;
        (*bufferC_st)[i] = Individual::characterLocus[i].C_st;
        (*bufferVarP)[i] = Individual::characterLocus[i].varP[0u];
        (*bufferVarG)[i] = Individual::characterLocus[i].varG[0u];
        (*bufferVarA)[i] = Individual::characterLocus[i].varA[0u];
        (*bufferVarD)[i] = Individual::characterLocus[i].varD;
        (*bufferVarI)[i] = Individual::characterLocus[i].varI[0u];
        (*bufferCovL)[i] = Individual::characterLocus[i].covL;
    }
    bufferFreq->flush();
    bufferF_it->flush();
    bufferF_is->flush();
    bufferF_st->flush();
    bufferP_st->flush();
    bufferG_st->flush();
    bufferQ_st->flush();
    bufferC_st->flush();
    bufferVarP->flush();
    bufferVarG->flush();
    bufferVarA->flush();
    bufferVarD->flush();
    bufferVarI->flush();
    bufferCovL->flush();

    // *** genome-wide decomposition of genetic variance (continued) ***
    // compute varA, varD, varI, covL and Fst by accumulating single locus contributions
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
        for(size_t cl = 0u; cl < 3u; ++cl) {
            Individual::varA[crctr][cl] = 0.0;
            Individual::varI[crctr][cl] = 0.0;
        }
        Individual::varD[crctr] = 0.0;
        Individual::covL[crctr] = 0.0;
        double var_t = 0.0, var_s = 0.0;
        for(size_t i : Individual::vertices[crctr]) {
            for(size_t cl = 0u; cl < 3u; ++cl) {
                Individual::varA[crctr][cl] += Individual::characterLocus[i].varA[cl];
                Individual::varI[crctr][cl] += Individual::characterLocus[i].varI[cl];
            }
            Individual::varD[crctr] += Individual::characterLocus[i].varD;
            Individual::covL[crctr] += Individual::characterLocus[i].covL;
            double p0 = Individual::characterLocus[i].alleleFrequency[0u];
            double p1 = Individual::characterLocus[i].alleleFrequency[1u];
            double p2 = Individual::characterLocus[i].alleleFrequency[2u];
            var_s += (n[1u] * sqr(p1) + n[2u] * sqr(p2)) / n[0u] - sqr(p0);
            var_t += p0 * (1.0 - p0);
        }
        Individual::F_st[crctr] = var_t < tiny ? 0.0 : var_s / var_t;
        Individual::P_st[crctr] = Xst(Individual::varP[crctr][0u],
                                      Individual::varP[crctr][1u],
                                      Individual::varP[crctr][2u], n);
        Individual::G_st[crctr] = Xst(Individual::varG[crctr][0u],
                                      Individual::varG[crctr][1u],
                                      Individual::varG[crctr][2u], n);
        Individual::Q_st[crctr] = Xst(Individual::varA[crctr][0u],
                                      Individual::varA[crctr][1u],
                                      Individual::varA[crctr][2u], n);
        Individual::C_st[crctr] =
            Xst(Individual::varD[crctr] + 0.5 * Individual::varI[crctr][0u],
                Individual::varI[crctr][1u],
                Individual::varI[crctr][2u], n);
    }
    recordData(t, n);
}

void analyseNetwork(int t)
{
    const char sep = ',';
    // *** node properties ***
    // open network link file
    std::ostringstream oss1;
    oss1 << "simulation_" << seed << '_' << t << "_nodes.csv";
    std::ofstream ofs(oss1.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    
    ofs << "id" << sep << "character" << sep << "linkage.group" << sep
        << "degree" << sep << "location" << sep
        << "effect.size" << sep << "dominance.coefficient" << sep
        << "allele.freq" << sep << "mean.effect" << sep
        << "avg.effect.substitution" << sep << "varP" << sep
        << "varG" << sep << "varA" << sep << "varD" << sep
        << "varI" << sep << "covL" << sep << "exp.heterozygosity" << sep
        << "Fit" << sep << "Fis" << sep << "Fst" << sep
        << "Gst" << sep << "Qst" << sep << "Cst" << '\n';
    
    for(size_t i = 0u; i < nLoci; ++i) {
        const size_t crctr = Individual::characterLocus[i].character;
        double pi = Individual::characterLocus[i].alleleFrequency[0u];
        ofs << i << sep
            << Individual::characterLocus[i].character << sep
            << Individual::characterLocus[i].linkageGroup << sep
            << Individual::characterLocus[i].edges.size() << sep
            << Individual::characterLocus[i].location << sep
            << scaleA[crctr] * Individual::characterLocus[i].effectSize << sep
            << scaleD[crctr] * Individual::characterLocus[i].dominanceCoeff << sep
            << pi << sep
            << Individual::characterLocus[i].meanEffect[0u] << sep
            << Individual::characterLocus[i].avgEffectOfSubstitution << sep
            << Individual::characterLocus[i].varP[0u] << sep
            << Individual::characterLocus[i].varG[0u] << sep
            << Individual::characterLocus[i].varA[0u] << sep
            << Individual::characterLocus[i].varD << sep
            << Individual::characterLocus[i].varI[0u] << sep
            << Individual::characterLocus[i].covL << sep
            << 2.0 * pi * (1.0 - pi) << sep
            << Individual::characterLocus[i].F_it << sep
            << Individual::characterLocus[i].F_is << sep
            << Individual::characterLocus[i].F_st << sep
            << Individual::characterLocus[i].G_st << sep
            << Individual::characterLocus[i].Q_st << sep
            << Individual::characterLocus[i].C_st <<'\n';
    }
    ofs.close();

    // *** edge properties ***
    // open network link file
    std::ostringstream oss2;
    oss2 << "simulation_" << seed << '_' << t << "_links.csv";
    ofs.open(oss2.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    ofs << "from" << sep << "to" << sep << "effect.size" << sep
        << "interaction.weight" << sep
        << "gene.freq.correlation" << sep
        << "genetic.correlation" << sep
        << "additive.genetic.correlation" << '\n';
    
    // determine edge properties
    const size_t n = population.size();
    for(size_t i = 0u; i < nLoci; ++i) {
        const size_t crctr = Individual::characterLocus[i].character;
        for(const std::pair<size_t, double> &edge : Individual::characterLocus[i].edges) {
            const size_t j = edge.first;
            const double eij = scaleI[crctr] * edge.second;
            double sumpipj = 0.0, sumgigj = 0.0, sumbibj = 0.0,
                sumxi = 0.0, sumxj = 0.0, sumxixi = 0.0, sumxjxj;
            double pi = Individual::characterLocus[i].alleleFrequency[0u];
            double pj = Individual::characterLocus[j].alleleFrequency[0u];
            double alphai = Individual::characterLocus[i].avgEffectOfSubstitution;
            double alphaj = Individual::characterLocus[j].avgEffectOfSubstitution;
            for(PInd pInd : population) {
                size_t ui = pInd->traitLocus[i].alleleCount;
                size_t uj = pInd->traitLocus[j].alleleCount;
                double xi = pInd->traitLocus[i].expression;
                double xj = pInd->traitLocus[j].expression;
                double gi = pInd->traitLocus[i].geneticValue;
                double gj = pInd->traitLocus[j].geneticValue;
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
            
            // variation in average effect dus to epistasis
            sumxi /= n;
            sumxj /= n;
            sumxixi = (sumxixi - n * sumxi * sumxi) / (n - 1u);
            sumxjxj = (sumxjxj - n * sumxj * sumxj) / (n - 1u);
            double ai = scaleA[crctr] * Individual::characterLocus[i].effectSize;
            double aj = scaleA[crctr] * Individual::characterLocus[j].effectSize;
            double weightij = fabs(aj) < tiny ? 0.0 :
                sqrt(sqr(Individual::characterLocus[j].avgEffectOfSubstitution * eij / aj)
                    * sumxixi);
            double weightji = fabs(ai) < tiny ? 0.0 :
                sqrt(sqr(Individual::characterLocus[i].avgEffectOfSubstitution * eij / ai)
                    * sumxjxj);
        
            // correlation between allele frequencies ~ population LD
            sumpipj = (sumpipj - n * pi * pj) / (n - 1u);
            double norm = sqrt(pi * (1.0 - pi) * pj * (1.0 - pj));
            double rp = norm > tiny ? sumpipj / norm : 0.0;
            
            // correlation between genetic values ~ genetic covariance
            double mui = Individual::characterLocus[i].meanEffect[0u];
            double muj = Individual::characterLocus[j].meanEffect[0u];
            double vari = Individual::characterLocus[i].varG[0u];
            double varj = Individual::characterLocus[j].varG[0u];
            sumgigj = (sumgigj - n * mui * muj) / (n - 1u);
            norm = sqrt(vari * varj);
            double rG = norm > tiny ? sumgigj / norm : 0.0;
            
            // correlation between breeding values ~ additive genetic covariance
            vari = Individual::characterLocus[i].varA[0u];
            varj = Individual::characterLocus[j].varA[0u];
            sumbibj /= (n - 1u);
            norm = sqrt(vari * varj);
            double rA = norm > tiny ? sumbibj / norm : 0.0;
            
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
    oss3 << "simulation_" << seed << '_' << t << "_dump.csv";
    ofs.open(oss3.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in networkAnalysis()");
    ofs << sep << "ecotype" << sep << "habitat" << sep << "attack.rate.1" << sep << "attack.rate.2";
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        ofs << sep << "trait." << crctr;
    ofs << '\n';
    size_t i = 0u;
    for(PInd pInd : population) {
        ofs << i  << sep << pInd->ecotype << sep << pInd->habitat
            << sep << pInd->attackRate.first << sep << pInd->attackRate.second;
        for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
            ofs << sep << pInd->traitP[crctr];
        ofs << '\n';
        ++i;
    }
    ofs.close();
}