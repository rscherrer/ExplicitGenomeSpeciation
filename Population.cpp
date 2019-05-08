//
// Created by p278834 on 7-5-2019.
//

#include "Population.h"
#include "random.h"
#include "queue"

void Population::dispersal(const ParameterSet& parameters)
{
    if(parameters.dispersalRate > 0.5) {
        for(PInd pInd : individuals)
            if(rnd::bernoulli(parameters.dispersalRate)) pInd -> disperse(parameters.nHabitat);
    }
    else {
        const size_t n = individuals.size();
        size_t k = rnd::binomial(n, parameters.dispersalRate);
        if(k == 0u) return;

        std::set<size_t> migrants;
        while(migrants.size() < k)
            migrants.insert(rnd::random_int(n));

        std::list<PInd>::const_iterator it = individuals.cbegin();
        size_t j = 0u;
        for(size_t i : migrants) {
            std::advance(it, i - j);
            (*it) -> disperse(parameters.nHabitat);
            j = i;
        }
    }
}


void Population::competitionAndReproduction(const size_t hab,
                                const ParameterSet& parameters,
                                std::vector<std::pair<double, double> >& resourceConsumption,
                                Individual::TradeOffPt& breakEvenPoint,
                                std::vector<std::pair<double, double> >& resourceEql,
                                std::vector<std::pair<size_t, size_t> >& genderCounts,
                                const Genome& genome,
                                const size_t nAccessibleResource = 2u)
{
    // accumulate attack rates and sort out females and males
    std::queue<PInd> females;
    std::vector<PInd> males;
    std::list<Individual::TradeOffPt> pts;
    std::list<PInd>::iterator iti = individuals.begin();
    resourceConsumption[hab].first = resourceConsumption[hab].second = 0.0;

    for(std::list<PInd>::iterator itj = individuals.end(); iti != itj;) {
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
    individuals.erase(individuals.begin(), iti);

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
                individuals.push_back(new Individual(fem, males[j], parameters, genome));
                if(parameters.costIncompat > 0.0) {
                    if(rnd::bernoulli(individuals.back()->getViability()))
                        individuals.pop_back();
                }
                --nOffspring;
            }
        }
        // female survival
        if(rnd::bernoulli(parameters.survivalProb)) individuals.push_back(fem);
        else delete fem;
    }
    // male survival
    for(size_t i = 0u; i < nm; ++i) {
        if(rnd::bernoulli(parameters.survivalProb)) individuals.push_back(males[i]);
        else delete males[i];
        males[i] = nullptr;
    }
}