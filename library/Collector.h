#ifndef EXPLICITGENOMESPECIATION_COLLECTOR_H
#define EXPLICITGENOMESPECIATION_COLLECTOR_H

// This header contains the declarations of the Collector class and allies
// This class is an analytical module for the simulation
// Using its main function "analyze", it computes many statistics from the
// population (genome-wide, locus-specific or edge-specific statistics)
// The content of the Collector is reset every time "analyze" is called

#include "Utilities.h"
#include "MetaPop.h"
#include "GenArch.h"
#include <cassert>

// A class for locus-specific statistics
class Locus
{

    friend class Collector;
    friend class Printer;
    friend class Connexion;

    Locus(const size_t &i, const size_t &t) :
        id(i),
        trait(t),
        gcounts(utl::uzeros(3u, 4u)),
        gsumgen(utl::zeros(3u, 4u)),
        gssqgen(utl::zeros(3u, 4u)),
        gbeta(std::vector<double>(3u, 0.0)),
        gexpec(std::vector<double>(3u, 0.0)),
        gmeans(std::vector<double>(3u, 0.0)),
        gdelta(std::vector<double>(3u, 0.0)),
        varG(std::vector<double>(3u, 0.0)),
        varP(std::vector<double>(3u, 0.0)),
        varA(std::vector<double>(3u, 0.0)),
        varN(std::vector<double>(3u, 0.0)),
        freqs(std::vector<double>(3u, 0.0)),
        varD(0.0),
        varI(0.0),
        varQ(0.0),
        varX(0.0),
        varS(0.0),
        Pst(0.0),
        Gst(0.0),
        Qst(0.0),
        Cst(0.0),
        Fst(0.0),
        alpha(0.0),
        meanQ(0.0),
        meanG(0.0),
        covQG(0.0),
        h(0.0),
        H(0.0),
        hobs(std::vector<double>(2u, 0.0)),
        tot(2u),
        all(3u),
        aa(0u),
        Aa(1u),
        AA(2u)
    {}

    void reset();
    void fillMatrices(const MetaPop&, const size_t&);
    void calcAlleleFreqs(const std::vector<size_t>&);
    void calcMeanQ();
    void calcMeanG(const size_t&);
    void calcVarQ(const size_t&);
    void calcCovQG(const size_t&);
    void calcAvgMutEffect();
    void calcVarX(const double&, const double&, const size_t&);
    void calcGenotypeStats();
    void calcVarG(const size_t&, const size_t&);
    void calcVarP(const size_t&, const double&);
    void calcVarA(const size_t&, const size_t&);
    void calcVarN(const size_t&, const size_t&);
    void calcVarD(const size_t&);
    void calcVarI(const size_t&);
    void calcVarS(const size_t&, const size_t&, const size_t&);
    void calcHWithin(const size_t&, const size_t&, const size_t&);
    void calcHAcross();
    void calcHObserved(const size_t&);
    void partitionVariance(const std::vector<size_t>&);

    size_t id; // locus number
    size_t trait; // encoded trait

    // Per ecotype per genotype
    std::vector<std::vector<size_t> > gcounts; // matrix of genotype counts
    std::vector<std::vector<double> > gsumgen; // matrix of sums of gen values
    std::vector<std::vector<double> > gssqgen; // matrix of ssqs of gen values

    // Note: the above matrices are used to compute all the subsequent
    // statistics in an efficient way and avoid unnecessary loops

    // per genotype
    std::vector<double> gbeta; // breeding values
    std::vector<double> gexpec; // expected additive genetic values
    std::vector<double> gmeans; // mean genetic values
    std::vector<double> gdelta; // dominance deviations

    // per ecotype
    std::vector<double> varG; // genetic variance
    std::vector<double> varP; // phenotypic variance
    std::vector<double> varA; // additive variance
    std::vector<double> varN; // non-additive variance
    std::vector<double> freqs; // allele frequencies

    double varD;
    double varI;
    double varQ;
    double varX;
    double varS;
    double Pst;
    double Gst;
    double Qst;
    double Cst;
    double Fst;
    double alpha;
    double meanQ;
    double meanG;
    double covQG;
    double h;
    double H;
    std::vector<double> hobs;

    // Indices for readability
    const size_t tot; // across ecotypes
    const size_t all; // across genotypes
    const size_t aa;
    const size_t Aa;
    const size_t AA;

};

// A class for edge-specific statistics
class Connexion
{
    friend class Collector;
    friend class Printer;

    Connexion(const size_t &e, const size_t &l0, const size_t &l1,
     const size_t &t) :
        id(e),
        i(l0),
        j(l1),
        trait(t),
        ggcounts(utl::uzeros(3u, 3u)),
        sprgen(0.0),
        corgen(0.0),
        corbreed(0.0),
        corfreq(0.0),
        avgi(0.0),
        avgj(0.0),
        tot(2u),
        all(3u),
        aa(0u),
        Aa(1u),
        AA(2u),
        bb(0u),
        Bb(1u),
        BB(2u)
    {}

    void reset();
    void fillMatrices(const MetaPop&, const size_t&);
    void calcCorGen(const size_t&, const std::vector<Locus>&);
    void calcCorBreed(const size_t&, const std::vector<Locus>&);
    void calcCorFreq(const size_t&, const std::vector<Locus>&);
    void calcBkgdIJ(const std::vector<Locus>&, const GenArch&, const Param&);

    double calcBkgdEffect(const Locus&, const Locus&, const double&,
     const double&, const double&, const double&) const;

    size_t id;
    size_t i;
    size_t j;
    size_t trait;

    std::vector<std::vector<size_t> > ggcounts;

    // Correlations in genetic values, breeding values and allele freq
    double sprgen;
    double corgen;
    double corbreed;
    double corfreq;

    // Variation in average effect due to epistasis
    double avgi;
    double avgj;

    // For readability
    const size_t tot;
    const size_t all; // across genotypes
    const size_t aa;
    const size_t Aa;
    const size_t AA;
    const size_t bb;
    const size_t Bb;
    const size_t BB;

};

typedef std::shared_ptr<std::ofstream> Stream;

// The collector class
class Collector
{

    friend class Printer;

public:

    Collector(const GenArch &arch) :
        counts(utl::uzeros(3u, 3u)),
        ecounts(std::vector<size_t>(3u, 0u)),
        sumgen(utl::zeros(3u, 3u, 3u)),
        sumphe(utl::zeros(3u, 3u, 3u)),
        ssqgen(utl::zeros(3u, 3u, 3u)),
        ssqphe(utl::zeros(3u, 3u, 3u)),
        means(utl::zeros(3u, 3u, 3u)),
        esumgen(utl::zeros(3u, 3u)),
        esumphe(utl::zeros(3u, 3u)),
        essqgen(utl::zeros(3u, 3u)),
        essqphe(utl::zeros(3u, 3u)),
        varG(utl::zeros(3u, 3u)),
        varP(utl::zeros(3u, 3u)),
        varA(utl::zeros(3u, 3u)),
        varN(utl::zeros(3u, 3u)),
        varD(std::vector<double>(3u, 0.0)),
        varI(std::vector<double>(3u, 0.0)),
        varT(std::vector<double>(3u, 0.0)),
        varS(std::vector<double>(3u, 0.0)),
        Pst(std::vector<double>(3u, 0.0)),
        Gst(std::vector<double>(3u, 0.0)),
        Qst(std::vector<double>(3u, 0.0)),
        Cst(std::vector<double>(3u, 0.0)),
        Fst(std::vector<double>(3u, 0.0)),
        genomescan(emptyloci(arch)),
        networkscan(emptyconnexions(arch)),
        EI(0.0),
        SI(0.0),
        RI(0.0),
        tot(2u)
    {}

    // Main function
    void analyze(const MetaPop&, const Param&, const GenArch&);

    // Various getters called in tests
    double getEI() const;
    double getSI() const;
    double getRI() const;
    double getVarP(const size_t&) const;
    double getVarN(const size_t&) const;
    double getVarI(const size_t&) const;
    double calcLocusGenotypeVarG(const size_t&, const size_t&, const MetaPop&,
     const size_t&) const;


    //----------------------------------
    // Functions by Thijs for the GUI
    // these getters are used in plotting
    // they are not optimized - they are not called that often.
    std::vector<double> get_Fst() const;
    std::vector<double> get_Gst() const;
    std::vector<double> get_eco_trait(const MetaPop &m) const;
    std::vector<double> get_eco_trait_deme(const MetaPop &m,
                                           size_t deme) const;
    std::vector<double> get_sex_trait(const MetaPop &m) const;
    std::vector<double> get_sex_trait_deme(const MetaPop &m,
                                           size_t deme) const;
    std::vector<double> get_neu_trait(const MetaPop &m) const;
    std::vector<double> get_neu_trait_deme(const MetaPop &m,
                                           size_t deme) const;

private:

    // Private makers
    std::vector<Locus> emptyloci(const GenArch&) const;
    std::vector<Connexion> emptyconnexions(const GenArch&) const;

    // Private setters
    void reset();
    void resetLocus(const size_t &l);
    void fillMatrices(const MetaPop&);
    void calcVarG(const size_t&, const size_t&);
    void calcVarP(const size_t&, const size_t&);
    void partitionVariance(const size_t&);
    void calcEI();
    void calcSI();
    void calcRI(const MetaPop&, const Param&);
    void analyzeLocus(const size_t&, const MetaPop&, const GenArch&, const Param&);
    void analyzeEdge(const size_t&, const MetaPop&, const GenArch&, const Param&);
    void analyzeTrait(const size_t&);

    // Private fields

    std::vector<std::vector<size_t> > counts; // per habitat per ecotype
    std::vector<std::size_t> ecounts; // per ecotype

    // per trait per habitat per ecotype
    std::vector<std::vector<std::vector<double> > > sumgen;
    std::vector<std::vector<std::vector<double> > > sumphe;
    std::vector<std::vector<std::vector<double> > > ssqgen;
    std::vector<std::vector<std::vector<double> > > ssqphe;
    std::vector<std::vector<std::vector<double> > > means;

    // per trait per ecotype
    std::vector<std::vector<double> > esumgen;
    std::vector<std::vector<double> > esumphe;
    std::vector<std::vector<double> > essqgen;
    std::vector<std::vector<double> > essqphe;

    // per trait per ecotype
    std::vector<std::vector<double> > varG;
    std::vector<std::vector<double> > varP;
    std::vector<std::vector<double> > varA;
    std::vector<std::vector<double> > varN;

    // per trait
    std::vector<double> varD;
    std::vector<double> varI;
    std::vector<double> varT;
    std::vector<double> varS;
    std::vector<double> Pst;
    std::vector<double> Gst;
    std::vector<double> Qst;
    std::vector<double> Cst;
    std::vector<double> Fst;

    std::vector<Locus> genomescan; // per locus
    std::vector<Connexion> networkscan; // per edge

    double EI;
    double SI;
    double RI;

    const size_t tot;

};

// Function to calculate F-statistics
double Xst(const std::vector<double>&, const std::vector<size_t>&);

#endif
