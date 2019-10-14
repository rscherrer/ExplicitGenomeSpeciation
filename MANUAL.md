# ExplicitGenomeSpeciation

An individual-based simulation of an adaptive speciation event, with explicit full genomes and quantitative genetics.

## Model

### Environment

There are two habitats connected by dispersal, labelled 0 and 1. Two resources are present across the lanscape, also labelled 0 and 1. We assume that resource are not necessarily equally distributed between the habitats, but can be asymmetrically distributed. Resource 0 is the most abundant resource in habitat 0, while resource 1 is most abundant in habitat 1. The concentration of the most abundant resource is the same in each habitat. The least abundant resource in each habitat has a concentration that is a proportion ```hsymmetry``` of that of the most abundant resource, where ```hsymmetry``` can range from 0 (resources are equally distributed between habitats) to 1 (there is only one resource available in each habitat). The concentration of the least abundant resource is also the same between the habitats.

### Life cycle

The dynamics of the population unfold in discrete time. At every time step, the following events occur in the corresponding order: dispersal, consumption, reproduction and survival.

#### Dispersal

Every individual can migrate to the alternative habitat with probability ```dispersal```. There is no dispersal during the burn-in period.

#### Consumption

We assume that the resource dynamics are fast relative to the population dynamics. We calculate the equilibrium concentrations of each resource in each habitat as a steady state of the resource dynamics. We provide two ways to model resource flows: logistic (```rdynamics = 0```) or chemostat (```rdynamics = 1```) dynamics. The logistic dynamics are characterized by the following differential equation:

![equation](img/logistic.jpg)

where *R* is the concentration of a given resource in a given habitat, *r* is the replenishment rate of the resource, *K* is the carrying capacity, and *C* is the consumption rate from the population. *C* is calculated as the sum of attack rates on this resource of all individuals in the present habitat. The chemostat dynamics are characterized by the equation:

![equation](img/chemostat.jpg)

where *a* is the rate of inflow of resource in the habitat and *b* the rate of outflow. *C* is calculated the same way as for the logistic dynamics.  

The fitness of each individual is calculated as:

![equation](img/fitness.jpg)

where *e<sub>0</sub>* and *e<sub>1</sub>* are the attack rates of the individual on each resource, given its ecological trait value *x*, and *R<sub>0</sub><sup>*</sup>* and *R<sub>1</sub><sup>*</sup>* are the equilibrium concentrations of the two resources in the individual's habitat. Note that only resource 0 is available during the burn-in period.

#### Reproduction

The fitness of individuals determines their reproductive success. At every generation, a mating seasons occurs for a number of discrete rounds, where females can choose males to mate with. It is possible that females remain unmated at the end of a mating season. The duration of the mating season is sampled each generation from a geometric distribution with parameter `matingcost`, which is the probability that the mating season will end after a given round. We recommend low values of `matingcost`, as to not strongly select against choosy females (who run high chances of not finding a mate if the mating season is too short).  

At each round of the mating season, females are presented with a male to assess. The probability of a male to encounter and being assessed by a females relative to other males is proportional to its fitness. Once encountered, the male is evaluated by the female, who can either accept or reject him as a mate. Once a female accepts a mate, the mating season ends for her and she produces a number of offspring sampled from a Poisson distribution with mean *f* `birth`, where `birth` is a baseline birth rate and *f* the fitness of the female.  

The probability to accept a given male depends on the female's mate preference *y* and on the male's ecological trait *x* value relative to her own. Positive values of *y* will favor assortative mating, i.e. with males that are more similar in ecological trait, while negative values will favor disassortative mating, i.e. with males that are more dissimilar in ecological trait. Mating is random if *y = 0*. Specifically, the mating probability of female *i* with male *j* is given, for *y > 0*, by:

![equation](img/assortative.jpg)

where... The mating probability in case *y < 0* is;

![equation](img/disassortative.jpg).

#### Survival

### Initialization

The founder population is populated with a certain number of individuals in each of the two habitats, defined by ```demesizes```. Founder individuals are generated with random genomes and therefore trait values may span a large range in the first time steps. This is not suitable for modelling e.g., a population that colonizes a new ecological niche, or resource, from a niche it is already adapted to. To model this scenario we incorporate a burn-in period in our simulation (of duration ```tburnin```), during which only resource 0 is available for use and there is no dispersal. If the population is initialized with all individuals in habitat 0, the burn-in period should canalize the population to specializing on resource 0, with ecological trait values gathered around -1. Once the burn-in period is over, the simulation runs for ```tend``` time steps where dispersal is allowed and both resources are available. If parameter ```record``` is set to 1, data are recorded every ```tsave``` time step outside the burn-in period.

### Life cycle

The simulation is run for ```tend``` discrete time steps.

## Parameters

Name | Definition
---|---
```rdynamics``` | whether the resources follow logistic (0) or a chemostat (1) dynamics
```inflow``` | rate of resource inflow in chemostat dynamics
```outflow``` | rate of resource outflow in chemostat dynamics
```capacity``` | resource carrying capacity in logistic dynamics
```replenish``` | resource growth rate in logistic dynamics
```hsymmetry``` | symmetry between the two habitats in resource carrying capacity or inflow rate
```ecosel``` | ecological selection coefficient
```dispersal``` | dispersal rate between habitats
```birth``` | expected number of offspring of a female with fitness 1
```survival``` | probability of an adult to survive to the next generation
```sexsel``` | strength of female mate preference
```matingcost``` | fitness cost of being choosy
```maxfeed``` | maximum achievable attack rate of an individual on a resource
```demesizes``` | initial population sizes of the two habitats
```nloci``` | number of loci in the genome
```nvertices``` | numbers of loci underlying each trait
```nedges``` | numbers of edges in the gene networks underlying each trait
```nchrom``` | number of chromosomes
```mutation``` | mutation rate per locus
```recombination``` | recombination rate (multiply by 100 to get genome size in centiMorgans)
```allfreq``` | frequency of allele 1 in the founder population
```scaleA``` | scaling parameter of additive effects for each trait
```scaleD``` | scaling parameter of dominance effects for each trait
```scaleI``` | scaling parameter of interaction effects for each trait
```scaleE``` | scaling parameter of environmental effects for each trait
```locusE``` | locus-specific environmental variance (initialize)
```skews``` | skew of the degree distributions of the gene networks underlying each trait
```effectshape``` | shape of the distribution of locus effect sizes
```effectscale``` | scale of the distribution of locus effect sizes
```interactionshape``` | shape of the distribution of gene interaction weights
```interactionscale``` | scale of the distribution of gene interaction weights
```dominancevar``` | variance of the distribution of locus dominance coefficients
```tburnin``` | duration of the burn-in period
```tend``` | duration of the simulation
```tsave``` | pace of data recording
```record``` | whether to record the data or not
```seed``` | seed of the random number generator
```ntrials``` | number of mating trials performed to compute mating isolation

