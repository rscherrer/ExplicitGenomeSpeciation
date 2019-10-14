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

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7BdR%7D%7Bdt%7D%20%3D%20r%20%5CBig%281%20-%20%5Cfrac%7BR%7D%7BK%7D%5CBig%29%20R%20%5C%2C%20-%20%5C%2C%20C%20%5C%2C%20R&bc=White&fc=Black&im=jpg&fs=12&ff=concmath&edit=0)

![equation](img/logistic.jpg)

where *R* is the concentration of a given resource in a given habitat, *r* is the replenishment rate of the resource, *K* is the carrying capacity, and *C* is the consumption rate from the population. *C* is calculated as the sum of attack rates on this resource of all individuals in the present habitat. The chemostat dynamics are characterized by the equation:

![eq2](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7BdR%7D%7Bdt%7D%20%3D%20a%20-%20b%20%5C%2C%20R%20-%20%5C%2C%20C%20%5C%2C%20R&bc=White&fc=Black&im=jpg&fs=12&ff=concmath&edit=0)

where *a* is the rate of inflow of resource in the habitat and *b* the rate of outflow. *C* is calculated the same way as for the logistic dynamics.  

The fitness of each individual is calculated as:

![img](http://www.sciweavers.org/tex2img.php?eq=f%20%3D%20e_0%28x%29%20%5C%2C%20R_0%5E%2A%20%2B%20e_1%28x%29%20%5C%2C%20R_1%5E%2A&bc=White&fc=Black&im=jpg&fs=12&ff=concmath&edit=0)

where *e_0* and *e_1* are the attack rates of the individual on each resource, given its ecological trait value *x*, and *R_0* and *R_1* are the equilibrium concentrations of the two resources in the individual's habitat. Note that only resource 0 is available during the burn-in period.

#### Reproduction

The fitness of individuals is taken as a measure of reproductive success.

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

