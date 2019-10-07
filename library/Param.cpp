#include "Param.h"


// Create a default seed based on clock
size_t Param::makeDefaultSeed()
{
    return static_cast<size_t>(std::chrono::high_resolution_clock::now().
     time_since_epoch().count());
}


// Make sure that the numbers of edges of the genetic networks do not go above
// their theoretical maximum, given the number of loci in each network
void Param::capEdges()
{
    for (size_t trait = 0u; trait < 3u; ++trait) {
        const size_t n = nvertices[trait];

        // Number of edges in a compete graph with N vertices
        const size_t emax = n * (n - 1u) / 2u;

        // Cap the number of edges
        if (nedges[trait] > emax) nedges[trait] = emax;
    }
}


// Functions to convert characters into integers
// Used to read parameters in

namespace fnv
{
  constexpr uint64_t _(uint64_t h, const char* s)
  {
    return (*s == 0) ? h :
      _((h * 1099511628211ull) ^ static_cast<uint64_t>(*s), s+1);
  }
}

constexpr uint64_t _(const char* s)
{
  return fnv::_(14695981039346656037ull, s);
}

uint64_t _(const std::string& s)
{
  return fnv::_(14695981039346656037ull, s.data());
}


// Read parameters from a file

void Param::read(const std::string &filename)
{
    std::ifstream inputfile;
    inputfile.open(filename);
    if (!inputfile.is_open()) {
        std::string msg = "Unable to open parameter file ";
        throw std::runtime_error(msg + filename);
    }

    update(inputfile);
    inputfile.close();
}

void Param::update(std::ifstream &file)
{

    std::string input;
    while (file >> input) {

        switch (_(input)) {

        case _("capacity"): file >> capacity; break;
        case _("replenish"): file >> replenish; break;
        case _("hsymmetry"): file >> hsymmetry; break;
        case _("ecosel"): file >> ecosel; break;
        case _("dispersal"): file >> dispersal; break;
        case _("birth"): file >> birth; break;
        case _("survival"): file >> survival; break;
        case _("sexsel"): file >> sexsel; break;
        case _("matingcost"): file >> matingcost; break;
        case _("maxfeed"): file >> maxfeed; break;
        case _("demesizes"):
            for (size_t i = 0u; i < 2u; ++i) file >> demesizes[i];
            break;
        case _("nvertices"):
            for (size_t i = 0u; i < 3u; ++i) file >> nvertices[i];
            break;
        case _("nedges"):
            for (size_t i = 0u; i < 3u; ++i) file >> nedges[i];
            break;
        case _("nchrom"): file >> nchrom; break;
        case _("mutation"): file >> mutation; break;
        case _("recombination"): file >> recombination; break;
        case _("allfreq"): file >> allfreq; break;
        case _("scaleA"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleA[i];
            break;
        case _("scaleD"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleD[i];
            break;
        case _("scaleI"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleI[i];
            break;
        case _("scaleE"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleE[i];
            break;
        case _("skews"):
            for (size_t i = 0u; i < 3u; ++i) file >> skews[i];
            break;
        case _("effectshape"): file >> effectshape; break;
        case _("effectscale"): file >> effectscale; break;
        case _("interactionshape"): file >> interactionshape; break;
        case _("interactionscale"): file >> interactionscale; break;
        case _("dominancevar"): file >> dominancevar; break;
        case _("tburnin"): file >> tburnin; break;
        case _("tend"): file >> tend; break;
        case _("tsave"): file >> tsave; break;
        case _("record"): file >> record; break;
        case _("seed"): file >> seed; break;

        default:
            throw std::runtime_error("Invalid parameter name: " + input); break;

        }
    }

    // Now update interactive parameters
    nloci = utl::sumu(nvertices);

    // Make sure genetic networks do not have more edges than feasible
    capEdges();

    // Check validity of parameter values
    checkParams();

    std::clog << "Parameters were read in succesfully.\n";

}


// Check that the parameter values are valid
void Param::checkParams()
{
    std::string msg = "No error detected";

    if (demesizes.size() != 2u)
        msg = "There should be two demes";
    if (dispersal < 0.0)
        msg = "Dispersal rate should be positive";
    if (dispersal > 1.0)
        msg = "Dispersal rate should be at most one";
    if (birth < 0.0)
        msg = "Birth rate should be positive";
    if (hsymmetry < 0.0)
        msg = "Habitat symmetry should be positive";
    if (hsymmetry > 1.0)
        msg = "Habitat symmetry should be at most one";
    if (survival < 0.0)
        msg = "Survival probability should be positive";
    if (survival > 1.0)
        msg = "Survival probability should be at most one";
    if (ecosel < 0.0)
        msg = "Selection coefficient should be positive";
    if (sexsel < 0.0)
        msg = "Mate preference strength should be positive";
    if (matingcost < 0.0)
        msg = "Mate evaluation cost should be positive";
    if (maxfeed < 0.0)
        msg = "Maximum feeding rate should be positive";
    if (capacity < 0.0)
        msg = "Maximum resource capacity should be positive";
    if (replenish < 0.0)
        msg = "Maximum resource growth should be positive";
    if (nvertices[0u] <= 1u)
        msg = "Numer of ecological loci should be at least two";
    if (nvertices[1u] <= 1u)
        msg = "Number of mating loci should be at least two";
    if (nvertices[2u] <= 1u)
        msg = "Number of neutral loci should be at least two";
    if (nchrom == 0u)
        msg = "Number of chromosomes should be at least one";
    if (nloci <= 5u)
        msg = "Total number of loci should be at least six";
    if (allfreq < 0.0)
        msg = "Frequency of SNPs should be positive";
    if (allfreq > 1.0)
        msg = "Frequency of SNPs should be at most one";
    if (mutation < 0.0)
        msg = "Mutation rate should be positive";
    if (mutation > 1.0)
        msg = "Mutation rate should be at most one";
    if (recombination < 0.0)
        msg = "Recombination rate should be positive";
    for (size_t i = 0u; i < 3u; ++i) {
        if (skews[i] < 0.0)
            msg = "Skewness should be positive";
        if (scaleA[i] < 0.0)
            msg = "Additive scaling should be positive";
        if (scaleD[i] < 0.0)
            msg = "Dominance scaling should be positive";
        if (scaleI[i] < 0.0)
            msg = "Interaction scaling should be positive";
        if (scaleE[i] < 0.0)
            msg = "Environmental scaling should be positive";
    }
    if (effectshape < 0.0)
        msg = "Effect size shape should be positive";
    if (effectscale < 0.0)
        msg = "Effect size scale should be positive";
    if (interactionshape < 0.0)
        msg = "Interaction weight shape should be positive";
    if (interactionscale < 0.0)
        msg = "Interaction weight scale should be positive";
    if (dominancevar < 0.0)
        msg = "Dominance variance should be positive";
    if (tburnin < 0)
        msg = "Burn-in time should be positive";
    if (tend <= 0)
        msg = "End time should be positive";
    if (tsave <= 0)
        msg = "Save time should be positive";

    if(msg != "No error detected")
        throw std::runtime_error(msg);
}
