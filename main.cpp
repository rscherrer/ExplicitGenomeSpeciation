/*==================================================================================================================================
                                                     main.cpp
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
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <cassert>
#include <gtest/gtest.h>
#include "random.h"
#include "individual.h"
#include "analysis.h"
#include "Genome.h"

/*=======================================================================================================
                                Parameter definitions and default values
========================================================================================================*/




/*=======================================================================================================
                                            I/O routines
========================================================================================================*/

template <class T>
bool read(const std::string &str, const std::string &name, T &par, std::ifstream &ifs)
{
    if(str == name) {
        ifs >> par;
        std::clog << "parameter " << str << " set to " << par << '\n';
        return true;
    }
    else return false;
}

void readParameters(const std::string& filename, ParameterSet& parameters)
{
    std::clog << "reading parameters from file " << filename << '\n';
    
    //open parameter file
    std::ifstream ifs(filename);
    if(!ifs.is_open())
        throw std::runtime_error("unable to open parameter file in readParameters()");
    
    std::string str;
    
    ifs >> str;
    if(str == "rng_seed_clock") parameters.seed = rnd::set_seed();
    else if(str == "rng_seed_user") {
        ifs >> parameters.seed;
        rnd::set_seed(parameters.seed);
    }
    else throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameterfile\n");
    
    ifs >> str;
    if(str == "architecture_generate") parameters.generateArchitecture = true;
    else if(str == "architecture_load") {
        parameters.generateArchitecture = false;
        ifs >> parameters.architecture;
    }
    else throw std::logic_error("\'architecture_generate\' or \'architecture_load <arg>\' expected at second line of parameterfile\n");
    
    while(ifs >> str) {
        if(read(str, "initial_sequence", parameters.sequenceString, ifs))
        {
            for(auto a : parameters.sequenceString)
            {
                parameters.sequence.push_back(a == '1');
            }
        }
        else if(str == "scale_A")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleA[crctr];
                std::clog   << "parameter " << "scale_A[" << crctr
                            << "] set to " << parameters.scaleA[crctr] << '\n';
            }
        else if(str == "scale_D")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleD[crctr];
                std::clog   << "parameter " << "scale_D[" << crctr
                            << "] set to " << parameters.scaleD[crctr] << '\n';
            }
        else if(str == "scale_I")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleI[crctr];
                std::clog   << "parameter " << "scale_I[" << crctr
                            << "] set to " << parameters.scaleI[crctr] << '\n';
            }
        else if(str == "scale_E")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleE[crctr];
                std::clog   << "parameter " << "scale_E[" << crctr
                            << "] set to " << parameters.scaleE[crctr] << '\n';
            }
        else if(read(str, "mutation_rate", parameters.mutationRate, ifs));
        else if(read(str, "genome_size_cm", parameters.mapLength, ifs));
        else if(read(str, "female_heterogamety", parameters.isFemaleHeteroGamety, ifs));
        else if(read(str, "dispersal_rate", parameters.dispersalRate, ifs));
        else if(read(str, "alpha", parameters.alpha, ifs));
        else if(read(str, "beta", parameters.beta, ifs));
        else if(read(str, "prob_survival", parameters.survivalProb, ifs));
        else if(read(str, "habitat_asymmetry", parameters.habitatAsymmetry, ifs));
        else if(read(str, "sel_coeff_ecol", parameters.ecoSelCoeff, ifs));
        else if(read(str, "preference_strength", parameters.matePreferenceStrength, ifs));
        else if(read(str, "preference_cost", parameters.mateEvaluationCost, ifs));
        else if(read(str, "incompatibility_cost", parameters.costIncompat, ifs));
        else if(read(str, "typeII_resource_utilisation", parameters.isTypeIIResourceUtilisation, ifs));
        else if(read(str, "typeII_mate_choice", parameters.isTypeIIMateChoice, ifs));
        else if(read(str, "network_skewness", parameters.networkSkewness, ifs));
        else if(str == "t_end") {
            ifs >> parameters.tBurnIn >> parameters.tEndSim;
            std::clog   << "burn-in period  " << parameters.tBurnIn << " generations \n";
            std::clog   << "simulation time " << parameters.tEndSim << " generations \n";
        }
        else if(str == "t_dat") {
            ifs >> parameters.tGetDat >> parameters.tSavDat;
            std::clog   << "data collected every " << parameters.tGetDat << " generations \n";
            std::clog   << "data stored    every " << parameters.tSavDat << " generations \n";
        }
        else throw std::runtime_error("unknown parameter " + str);
    }
    std::clog << "parameters were read in successfully\n";
}

void writeParameters(std::ofstream &ofs, const ParameterSet& parameters, const char sep = ' ')
{
    ofs << "rng_seed"       << sep  << parameters.seed         << '\n';
    ofs << "architecture"   << sep  << parameters.architecture << '\n';
    ofs << "scale_A";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleA[crctr];
    ofs << '\n';
    ofs << "scale_D";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleD[crctr];
    ofs << '\n';
    ofs << "scale_I";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleI[crctr];
    ofs << '\n';
    ofs << "scale_E";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleE[crctr];
    ofs << '\n';
    ofs << "mutation_rate"  << sep  << parameters.mutationRate << '\n';
    ofs << "genome_size_cm" << sep  << parameters.mapLength    << '\n';
    ofs << "female_heterogamety" << sep  << parameters.isFemaleHeteroGamety << '\n';
    ofs << "dispersal_rate" << sep  << parameters.dispersalRate << '\n';
    ofs << "alpha" << sep  << parameters.alpha << '\n';
    ofs << "beta" << sep  << parameters.beta << '\n';
    ofs << "prob_survival" << sep  << parameters.survivalProb << '\n';
    ofs << "habitat_asymmetry" << sep  << parameters.habitatAsymmetry << '\n';
    ofs << "sel_coeff_ecol" << sep  << parameters.ecoSelCoeff << '\n';
    ofs << "preference_strength" << sep  << parameters.matePreferenceStrength << '\n';
    ofs << "preference_cost" << sep  << parameters.mateEvaluationCost << '\n';
    ofs << "incompatibility_cost" << sep << parameters.costIncompat << '\n';
    ofs << "typeII_resource_utilisation" << sep << parameters.isTypeIIResourceUtilisation << '\n';
    ofs << "typeII_mate_choice" << sep << parameters.isTypeIIMateChoice << '\n';
    ofs << "network_skewness" << sep << parameters.networkSkewness << '\n';
    ofs << "t_end" << sep  << parameters.tBurnIn << sep << parameters.tEndSim << '\n';
    ofs << "t_dat" << sep  << parameters.tGetDat << sep << parameters.tSavDat << '\n';
    ofs << "initial_sequence" << sep << (parameters.sequence.size() == parameters.nBits ? to_string(parameters.sequence) : "random") << '\n';
}




/*=======================================================================================================
                                    biological model implementation
========================================================================================================*/

void dispersal(std::list<PInd>& population, const ParameterSet& parameters)
{
    if(parameters.dispersalRate > 0.5) {
        for(PInd pInd : population)
            if(rnd::bernoulli(parameters.dispersalRate)) pInd -> disperse(parameters.nHabitat);
    }
    else {
        const size_t n = population.size();
        size_t k = rnd::binomial(n, parameters.dispersalRate);
        if(k == 0u) return;
    
        std::set<size_t> migrants;
        while(migrants.size() < k)
            migrants.insert(rnd::random_int(n));
    
        std::list<PInd>::const_iterator it = population.cbegin();
        size_t j = 0u;
        for(size_t i : migrants) {
            std::advance(it, i - j);
            (*it) -> disperse(parameters.nHabitat);
            j = i;
        }
    }
}


bool tradeOffCompare (const Individual::TradeOffPt &x, const Individual::TradeOffPt &y) {
    bool yOnLeft = y.first < y.second;
    if(x.first < x.second) {
        if(yOnLeft) return (x.first < y.first);
        else return true;
    }
    else {
        if(yOnLeft) return false;
        else return (x.second < y.second);
    }
}

void competitionAndReproduction(const size_t hab,
                                const ParameterSet& parameters,
                                std::list<PInd>& population,
                                std::vector<std::pair<double, double> >& resourceConsumption,
                                Individual::TradeOffPt& breakEvenPoint,
                                std::vector<std::pair<double, double> >& resourceEql,
                                std::vector<std::pair<size_t, size_t> >& genderCounts,
                                const size_t nAccessibleResource = 2u)
{
    // accumulate attack rates and sort out females and males
    std::queue<PInd> females;
    std::vector<PInd> males;
    std::list<Individual::TradeOffPt> pts;
    std::list<PInd>::iterator iti = population.begin();
    resourceConsumption[hab].first = resourceConsumption[hab].second = 0.0;
        
    for(std::list<PInd>::iterator itj = population.end(); iti != itj;) {
        if((*iti)->getHabitat() == hab) {
            Individual::TradeOffPt pt = (*iti)->getAttackRate();
            if(nAccessibleResource < 2u) pt.second = 0.0;
            pts.push_back(pt);
            
            // for the moment, assume that the individual utilises the first resource
            resourceConsumption[hab].first += pt.first;

            // but sum the attack rates on the second resource anyway if type II resource utilisation
            if(parameters.isTypeIIResourceUtilisation) resourceConsumption[hab].second += pt.second;

            if((*iti)->isFemale(parameters.isFemaleHeteroGamety)) females.push(*iti);
            else males.push_back(*iti);
            ++iti;
        }
        else { // move individuals in the other habitat towards the end
            --itj;
            std::swap(*iti, *itj);
        }
    }
    population.erase(population.begin(), iti);
    
    // determine equilibrium scaled resource densities and assign final ecotype

    pts.sort(tradeOffCompare);
    breakEvenPoint = pts.back();
    breakEvenPoint.first *= 0.5;
    breakEvenPoint.second = 0.0;

    // security check
    if(!parameters.isTypeIIResourceUtilisation) if(resourceConsumption[hab].second != 0.0) throw std::logic_error("consumption of the second resource should be zero");

    // find resource equilibrium and break-even point (used only for ecotype classification in type II resource utilisation)
    resourceEql[hab].first = (hab == 0u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].first);
    resourceEql[hab].second = (hab == 1u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].second);

    for(const Individual::TradeOffPt &pt : pts) {
        if(!parameters.isTypeIIResourceUtilisation) {
            resourceEql[hab].first = (hab == 0u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].first);
            resourceEql[hab].second = (hab == 1u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].second);
        }
        if(pt.first * resourceEql[hab].first < pt.second * resourceEql[hab].second) {
            if(!parameters.isTypeIIResourceUtilisation) {
                // switching from resource 1 to 2 is beneficial
                resourceConsumption[hab].first -= pt.first;
                resourceConsumption[hab].second += pt.second;
            }
        }
        else {
            // set break-even point to be used in later ecotype classification
            breakEvenPoint = pt;
            break;
        }
    }

    //std::cout << hab << " : " << resourceEql[hab].first << ' ' << resourceEql[hab].second << '\n';
    
    // mate choice and offspring production
    const size_t nf = genderCounts[hab].first = females.size();
    const size_t nm = genderCounts[hab].second = males.size();
    
    // terminate if there are no males or no females
    if(nf == 0u || nm == 0u) return;
    
    // compute reproductive success for males
    std::vector<double> maleSuccess(nm);
    double sum = 0.0;
    for(size_t i = 0u; i < nm; ++i) {
        // pick the resource that yields the highest payoff (not if type II resource utilisation)
        Individual::TradeOffPt pt = males[i]->getAttackRate();
        if(nAccessibleResource < 2u) pt.second = 0.0;
        maleSuccess[i] = parameters.isTypeIIResourceUtilisation ? pt.first * resourceEql[hab].first + pt.second * resourceEql[hab].second : std::max(pt.first * resourceEql[hab].first, pt.second * resourceEql[hab].second);
        // add stabilising selection on mating trait during burn-in period
        if(nAccessibleResource < 2u)
            maleSuccess[i] *= males[i]->getBurnInRpSc(parameters.ecoSelCoeff);
        sum += maleSuccess[i];
        
    }
    if(sum < parameters.tiny) maleSuccess = std::vector<double>(nm, 1.0);
    std::discrete_distribution<size_t> maleMarket(maleSuccess.begin(), maleSuccess.end());

    // sample family sizes for females and implement mate choice
    const size_t seasonEnd = rnd::geometric(parameters.mateEvaluationCost);
    while(!females.empty())
    {
        PInd fem = females.front();
        females.pop();
    
        // compute female reproductive success
        Individual::TradeOffPt pt = fem->getAttackRate();
        if(nAccessibleResource < 2u) pt.second = 0.0;
        double femaleSuccess = parameters.isTypeIIResourceUtilisation ? pt.first * resourceEql[hab].first + pt.second * resourceEql[hab].second : std::max(pt.first * resourceEql[hab].first, pt.second * resourceEql[hab].second);
        if(nAccessibleResource < 2u)
            femaleSuccess *= fem->getBurnInRpSc(parameters.ecoSelCoeff);
    
        // sample family size for female
        size_t nOffspring = rnd::poisson(parameters.beta * femaleSuccess);
        fem->prepareChoice();
        for(size_t t = 0u; nOffspring && t < seasonEnd; ++t) {
            
            // sample a male
            const size_t j = maleMarket(rnd::rng);

            if(fem->acceptMate(males[j], parameters)) {
                // add offspring to the population only if it survives development
                population.push_back(new Individual(fem, males[j], parameters));
                if(parameters.costIncompat > 0.0) {
                    if(rnd::bernoulli(population.back()->getViability()))
                        population.pop_back();
                }
                --nOffspring;
            }
        }
        // female survival
        if(rnd::bernoulli(parameters.survivalProb)) population.push_back(fem);
        else delete fem;
    }
    // male survival
    for(size_t i = 0u; i < nm; ++i) {
        if(rnd::bernoulli(parameters.survivalProb)) population.push_back(males[i]);
        else delete males[i];
        males[i] = nullptr;
    }
}




/*=======================================================================================================
                                            main()
========================================================================================================*/

int main(int argc, char * argv[])
{

    #if !defined(NDEBUG)

        // In debug mode run tests only
        std::cout << "Running tests in debug mode...";
        testing::InitGoogleTest(&argc, argv);
        int testResult = RUN_ALL_TESTS();

    #else

    std::ofstream logFile; // must be accessible from catch

    try {

        // *** preliminaries ***

        // Parameters
        ParameterSet parameters;

        std::ofstream datFile, arcFile;

        // Declare pointers to addresses where buffers are stored
        // The values (the buffers) are accessed with (*pointerName)
        Buffer *bufferFreq, *bufferF_it, *bufferF_is, *bufferF_st,
                *bufferP_st, *bufferG_st, *bufferQ_st, *bufferC_st,
                *bufferVarP, *bufferVarG, *bufferVarA, *bufferVarD, *bufferVarI;

        // set parameter values
        if(argc == 1) {
            parameters.seed = rnd::set_seed(); // use default parameters values and use clock to set random seed
            parameters.generateArchitecture = true;
        }
        else if(argc == 2)
            readParameters(argv[1], parameters);
        else throw std::runtime_error("invalid number of program arguments in main()");

        Population population;
        Genome genome;

        // initialise genetic architecture
        if(parameters.generateArchitecture) {
            std::ostringstream oss;
            oss << "architecture_" << parameters.seed << ".txt";
            parameters.architecture = oss.str();
            genome.generateGeneticArchitecture(parameters);
            genome.storeGeneticArchitecture(parameters.architecture, parameters);
        }
        else genome.loadGeneticArchitecture(parameters.architecture, parameters, population);



        // open files and data buffers
        std::clog << "opening files and data buffers.";
        std::ostringstream oss;
        oss << "simulation_" << parameters.seed;
        logFile.open(oss.str() + ".log");
        datFile.open(oss.str() + ".dat");
        arcFile.open(oss.str() + "_fossil_record.txt");
        if(!(logFile.is_open() && datFile.is_open() && arcFile.is_open()))
            throw std::runtime_error("unable to open output files in main()");
        BufferBox bufferPointers;
        bufferPointers.bufferFreq = new Buffer("freq", parameters);
        bufferPointers.bufferF_it = new Buffer("Fit", parameters);
        bufferPointers.bufferF_is = new Buffer("Fis", parameters);
        bufferPointers.bufferF_st = new Buffer("Fst", parameters);
        bufferPointers.bufferP_st = new Buffer("Pst", parameters);
        bufferPointers.bufferG_st = new Buffer("Gst", parameters);
        bufferPointers.bufferQ_st = new Buffer("Qst", parameters);
        bufferPointers.bufferC_st = new Buffer("Cst", parameters);
        bufferPointers.bufferVarP = new Buffer("varP", parameters);
        bufferPointers.bufferVarG = new Buffer("varG", parameters);
        bufferPointers.bufferVarA = new Buffer("varA", parameters);
        bufferPointers.bufferVarD = new Buffer("varD", parameters);
        bufferPointers.bufferVarI = new Buffer("varI", parameters);
        std::clog << "..done\n";

        // store parameter values
        std::clog << "storing parameter values..";
        logFile << "parameters: ";
        if(argc == 1) logFile << "default values\n";
        else logFile << "imported from file " << argv[1] << '\n';
        writeParameters(logFile, parameters);
        std::clog << "..done\n";

        // write data file header
        std::clog << "writing data file header.";
        datFile << '\t' << "pop.size"
                << '\t' << "females"
                << '\t' << "males";
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "pop.size." << hab;
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "attack.rate." << hab << ".1" << '\t' << "attack.rate." << hab << ".2";
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "resource." << hab << ".1" << '\t' << "resource." << hab << ".2";
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
            datFile << '\t' << "pop.1.size" << crctr
                    << '\t' << "pop.2.size" << crctr
                    << '\t' << "mean.all." << crctr
                    << '\t' << "mean.1." << crctr
                    << '\t' << "mean.2." << crctr
                    << '\t' << "varP." << crctr
                    << '\t' << "varG." << crctr
                    << '\t' << "varA." << crctr
                    << '\t' << "varD." << crctr
                    << '\t' << "varI." << crctr
                    << '\t' << "Fst."  << crctr
                    << '\t' << "Pst."  << crctr
                    << '\t' << "Gst."  << crctr
                    << '\t' << "Qst."  << crctr
                    << '\t' << "Cst."  << crctr;
        datFile << '\t' << "speciation.cube.spatial.isolation"
                << '\t' << "speciation.cube.ecological.isolation"
                << '\t' << "speciation.cube.mating.isolation"
                << '\t' << "post.zygotic.isolation"<< '\n';
        std::clog << "..done\n";





        // record start of simulation
        auto tStart = std::chrono::system_clock::now();

        // *** simulation ***
        // create initial population

        std::list<PInd> population;
        std::vector<std::pair<double, double> > resourceConsumption, resourceEql;
        std::vector<std::pair<size_t, size_t> > genderCounts;
        Individual::TradeOffPt breakEvenPoint;

        std::clog << "creating initial population.";
        if(parameters.sequence.size() == parameters.nBits)
            for(size_t i = 0u; i < parameters.nIndividualInit ; ++i)
                population.push_back(new Individual(parameters.sequence, parameters));
        else for(size_t i = 0u; i < parameters.nIndividualInit ; ++i)
                population.push_back(new Individual(parameters));
        std::clog << "..done\n";

        // enter simulation loop
        std::clog << "entering simulation loop:\n";
        for(int t = 1 - parameters.tBurnIn; t <= parameters.tEndSim; ++t) {
            if(t > 0) {
                // default
                dispersal(population, parameters);
                competitionAndReproduction(0u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts);
                competitionAndReproduction(1u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts);

            }
            else competitionAndReproduction(0u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts, 1u); // burn-in period
            if(population.size() < 2u) {
                std::clog << "population size underflow at t = " << t << '\n';
                break;
            }
            if(t % parameters.tGetDat == 0u) decomposeVariance(t, parameters, bufferPointers, breakEvenPoint, resourceConsumption, resourceEql, arcFile, datFile, genderCounts, population, genome);
            if(t % parameters.tSavDat == 0u) analyseNetwork(t, parameters, population, genome);
        }

        // *** finalisation ***
        // record end of simulation
        auto tEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = tEnd-tStart;
        logFile << "Time to complete simulation: " << diff.count() << " s\n";

        // close output file
        logFile.close();
        datFile.close();
        arcFile.close();

        // free allocated memory
        while(!population.empty()) {
            delete population.back();
            population.pop_back();
        }
    }
    catch(const std::exception &err) {
        std::cerr << "exception: " << err.what() << '\n';
        logFile << "exception: " << err.what() << '\n';
        exit(EXIT_FAILURE);
    }

    #endif


    return EXIT_SUCCESS;
}