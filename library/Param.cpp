#include "Param.h"

// Constructor
//------------

Param::Param() :
    rdynamics(1u),
    capacity(1.0),
    replenish(2375.0),
    inflow(400.0),
    outflow(100.0),
    hsymmetry(0.0),
    ecosel(1.8),
    dispersal(1.0E-2),
    birth(1.0),
    survival(0.8),
    sexsel(10.0),
    matingcost(0.01),
    ecoscale(1.0),
    demesizes({ 100u, 0u }),
    nloci(90u), // cannot be provided
    nvertices({ 30u, 30u, 30u }),
    nedges({ 30u, 30u, 30u }),
    nchrom(3u),
    mutation(1.0E-3),
    recombination(3.0),
    allfreq(0.2),
    scaleA({ 1.0, 1.0, 1.0 }),
    scaleD({ 0.0, 0.0, 0.0 }),
    scaleI({ 0.0, 0.0, 0.0 }),
    scaleE({ 0.0, 0.0, 0.0 }),
    locusE({ 0.0, 0.0, 0.0 }), // cannot be provided
    skews({ 1.0, 1.0, 1.0 }),
    effectshape(2.0),
    effectscale(1.0),
    interactionshape(5.0),
    interactionscale(1.0),
    dominancevar(1.0),
    tburnin(0),
    tend(10),
    tsave(10),
    tfreeze(100),
    tpedigree(15000),
    talkative(true),
    record(true),
    datsave(true),
    choosewhattosave(false),
    gensave(false),
    archsave(false),
    archload(false),
    parsave(true),
    pedigreesave(false),
    archfile("architecture.txt"),
    parfile("paramlog.txt"),
    orderfile("whattosave.txt"),
    logfile("log.txt"),
    freezerfile("freezer.dat"),
    locifile("locivalues.dat"),
    pedigreefile("pedigree.dat"),
    seed(makeDefaultSeed()),
    ntrials(100u),
    pedigreetrials(100u),
    pedigreeoffspring(10u)
{

    // Make sure parameter values make sense
    check();

    // Seed the random number generator
    rnd::rng.seed(seed);
}

// Member functions
//-----------------

// Create a default seed based on clock
size_t Param::makeDefaultSeed()
{
    return static_cast<size_t>(std::chrono::high_resolution_clock::now().
     time_since_epoch().count());
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

    import(inputfile);
    inputfile.close();
}

void Param::import(std::ifstream &file)
{

    std::string input;

    while (file >> input) {

        if (input == "rdynamics") file >> rdynamics;
        else if (input == "capacity") file >> capacity;
        else if (input == "replenish") file >> replenish;
        else if (input == "inflow") file >> inflow;
        else if (input == "outflow") file >> outflow;
        else if (input == "hsymmetry") file >> hsymmetry;
        else if (input == "ecosel") file >> ecosel;
        else if (input == "dispersal") file >> dispersal;
        else if (input == "birth") file >> birth;
        else if (input == "survival") file >> survival;
        else if (input == "sexsel") file >> sexsel;
        else if (input == "matingcost") file >> matingcost;
        else if (input == "ecoscale") file >> ecoscale;
        else if (input == "demesizes")
            for (size_t i = 0u; i < 2u; ++i) file >> demesizes[i];
        else if (input == "nvertices")
            for (size_t i = 0u; i < 3u; ++i) file >> nvertices[i];
        else if (input == "nedges")
            for (size_t i = 0u; i < 3u; ++i) file >> nedges[i];
        else if (input == "nchrom") file >> nchrom;
        else if (input == "mutation") file >> mutation;
        else if (input == "recombination") file >> recombination;
        else if (input == "allfreq") file >> allfreq;
        else if (input == "scaleA")
            for (size_t i = 0u; i < 3u; ++i) file >> scaleA[i];
        else if (input == "scaleD")
            for (size_t i = 0u; i < 3u; ++i) file >> scaleD[i];
        else if (input == "scaleI")
            for (size_t i = 0u; i < 3u; ++i) file >> scaleI[i];
        else if (input == "scaleE")
            for (size_t i = 0u; i < 3u; ++i) file >> scaleE[i];
        else if (input == "skews")
            for (size_t i = 0u; i < 3u; ++i) file >> skews[i];
        else if (input == "effectshape") file >> effectshape;
        else if (input == "effectscale") file >> effectscale;
        else if (input == "interactionshape") file >> interactionshape;
        else if (input == "interactionscale") file >> interactionscale;
        else if (input == "dominancevar") file >> dominancevar;
        else if (input == "tburnin") file >> tburnin;
        else if (input == "tend") file >> tend;
        else if (input == "tsave") file >> tsave;
        else if (input == "tfreeze") file >> tfreeze;
        else if (input == "tpedigree") file >> tpedigree;
        else if (input == "talkative") file >> talkative;
        else if (input == "record") file >> record;
        else if (input == "datsave") file >> datsave;
        else if (input == "choosewhattosave") file >> choosewhattosave;
        else if (input == "gensave") file >> gensave;
        else if (input == "archsave") file >> archsave;
        else if (input == "archload") file >> archload;
        else if (input == "parsave") file >> parsave;
        else if (input == "pedigreesave") file >> pedigreesave;
        else if (input == "archfile") file >> archfile;
        else if (input == "parfile") file >> parfile;
        else if (input == "orderfile") file >> orderfile;
        else if (input == "logfile") file >> logfile;
        else if (input == "freezerfile") file >> freezerfile;
        else if (input == "locifile") file >> locifile;
        else if (input == "pedigreefile") file >> pedigreefile;
        else if (input == "seed") file >> seed;
        else if (input == "ntrials") file >> ntrials;
        else if (input == "pedigreetrials") file >> pedigreetrials;
        else if (input == "pedigreeoffspring") file >> pedigreeoffspring;
        else
            throw std::runtime_error("Invalid parameter name: " + input);

    }

    // Now update interactive parameters
    update();

    std::clog << "Parameters were read in succesfully.\n";

}

void Param::update()
{
    rnd::rng.seed(seed);
    nloci = utl::sum(nvertices);
    check();
}

// Check that the parameter values are valid
void Param::check() const
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
    if (ecoscale < 0.0)
        msg = "Ecological scale should be positive";
    if (capacity <= 0.0)
        msg = "Maximum resource capacity should be positive";
    if (replenish <= 0.0)
        msg = "Maximum resource growth should be positive";
    if (rdynamics > 1u)
        msg = "Resource dynamics is either 0 (logistic) or 1 (chemostat)";
    if (inflow <= 0.0)
        msg = "Resource inflow rate should be positive";
    if (outflow <= 0.0)
        msg = "Resource outflow rate should be positive";
    for (size_t i = 0u; i < 3u; ++i)
        if (nvertices[i] < 2u)
            msg = "Number of loci per trait should be at least two";
    for (size_t i = 0u; i < 3u; ++i) {
        const size_t n = nvertices[i];
        const bool cond = nedges[i] >= n - 1u && nedges[i] <= n * (n - 1u) / 2u;
        if (!cond)
            msg = "Number of edges per trait should be between n-1 and n(n-1)/2";
    }
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
    if (tfreeze <= 0)
        msg = "Freezing time should be positive";
    if (ntrials == 0u)
        msg = "Number of mating trials should be at least one";

    if (msg != "No error detected")
        throw std::runtime_error(msg);
}

void Param::save() const
{
    std::ofstream file(parfile);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file " + parfile);
    write(file);
    file.close();
}

void Param::write(std::ofstream &file) const
{

    file << "rdynamics " << rdynamics << '\n';
    file << "capacity " << capacity << '\n';
    file << "replenish " << replenish << '\n';
    file << "inflow " << inflow << '\n';
    file << "outflow " << outflow << '\n';
    file << "hsymmetry " << hsymmetry << '\n';
    file << "ecosel " << ecosel << '\n';
    file << "dispersal " << dispersal << '\n';
    file << "birth " << birth << '\n';
    file << "survival " << survival << '\n';
    file << "sexsel " << sexsel << '\n';
    file << "matingcost " << matingcost << '\n';
    file << "ecoscale " << ecoscale << '\n';
    file << "demesizes ";
    for (size_t i = 0u; i < 2u; ++i) file << demesizes[i] << ' ';
    file << '\n';
    file << "nvertices ";
    for (size_t i = 0u; i < 3u; ++i) file << nvertices[i] << ' ';
    file << '\n';
    file << "nedges ";
    for (size_t i = 0u; i < 3u; ++i) file << nedges[i] << ' ';
    file << '\n';
    file << "nchrom " << nchrom << '\n';
    file << "mutation " << mutation << '\n';
    file << "recombination " << recombination << '\n';
    file << "allfreq " << allfreq << '\n';
    file << "scaleA ";
    for (size_t i = 0u; i < 3u; ++i) file << scaleA[i] << ' ';
    file << '\n';
    file << "scaleD ";
    for (size_t i = 0u; i < 3u; ++i) file << scaleD[i] << ' ';
    file << '\n';
    file << "scaleI ";
    for (size_t i = 0u; i < 3u; ++i) file << scaleI[i] << ' ';
    file << '\n';
    file << "scaleE ";
    for (size_t i = 0u; i < 3u; ++i) file << scaleE[i] << ' ';
    file << '\n';
    file << "skews ";
    for (size_t i = 0u; i < 3u; ++i) file << skews[i] << ' ';
    file << '\n';
    file << "effectshape " << effectshape << '\n';
    file << "effectscale " << effectscale << '\n';
    file << "interactionshape " << interactionshape << '\n';
    file << "interactionscale " << interactionscale << '\n';
    file << "dominancevar " << dominancevar << '\n';
    file << "tburnin " << tburnin << '\n';
    file << "tend " << tend << '\n';
    file << "tsave " << tsave << '\n';
    file << "tfreeze " << tfreeze << '\n';
    file << "tpedigree " << tpedigree << '\n';
    file << "talkative " << talkative << '\n';
    file << "record " << record << '\n';
    file << "datsave " << datsave << '\n';
    file << "choosewhattosave " << choosewhattosave << '\n';
    file << "gensave " << gensave << '\n';
    file << "archsave " << archsave << '\n';
    file << "archload " << archload << '\n';
    file << "parsave " << parsave << '\n';
    file << "pedigreesave " << pedigreesave << '\n';
    file << "archfile " << archfile << '\n';
    file << "parfile " << parfile << '\n';
    file << "orderfile " << orderfile << '\n';
    file << "logfile " << logfile << '\n';
    file << "freezerfile " << freezerfile << '\n';
    file << "locifile " << locifile << '\n';
    file << "pedigreefile " << pedigreefile << '\n';
    file << "seed " << seed << '\n';
    file << "ntrials " << ntrials << '\n';
    file << "pedigreetrials " << pedigreetrials << '\n';
    file << "pedigreeoffspring " << pedigreeoffspring << '\n';

}
