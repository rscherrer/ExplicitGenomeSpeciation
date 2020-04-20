#ifndef EXPLICITGENOMESPECIATION_COLLECTOR_H
#define EXPLICITGENOMESPECIATION_COLLECTOR_H


#include "Utilities.h"
#include "MetaPop.h"
#include "GenArch.h"
#include <cassert>

struct Locus
{
    Locus(const size_t &i, const size_t &t) :
        id(i),
        trait(t),
        varG(utl::zeros(3u)),
        varP(utl::zeros(3u)),
        varA(utl::zeros(3u)),
        varN(utl::zeros(3u)),
        varD(0.0),
        varI(0.0),
        varZ(0.0),
        varX(0.0),
        Pst(0.0),
        Gst(0.0),
        Qst(0.0),
        Cst(0.0),
        Fst(0.0),
        alpha(0.0),
        beta(utl::zeros(3u)),
        meang(0.0),
        freq(0.0)
    {}

    size_t id;
    size_t trait;

    std::vector<double> varG; // per ecotype
    std::vector<double> varP; // per ecotype
    std::vector<double> varA; // per ecotype
    std::vector<double> varN; // per ecotype
    double varD;
    double varI;
    double varZ;
    double varX;

    double Pst;
    double Gst;
    double Qst;
    double Cst;
    double Fst;

    double alpha;
    std::vector<double> beta; // per genotype
    double meang;
    double freq;
};

struct Connexion
{
    Connexion(const size_t &e, const size_t &i, const size_t &j,
     const size_t &t) :
        id(e),
        loc1(i),
        loc2(j),
        trait(t),
        corgen(0.0),
        corbreed(0.0),
        corfreq(0.0),
        avgi(0.0),
        avgj(0.0)
    {}

    size_t id;
    size_t loc1;
    size_t loc2;
    size_t trait;

    // Correlations in genetic values, breeding values and allele freq
    double corgen;
    double corbreed;
    double corfreq;

    // Variation in average effect due to epistasis
    double avgi;
    double avgj;
};

typedef std::shared_ptr<std::ofstream> Stream;

class Collector
{

    friend class Printer;

public:

    Collector(const GenArch &arch) :
        counts(utl::uzeros(3u, 3u)),
        means(utl::zeros(3u, 3u, 3u)),
        varG(utl::zeros(3u, 3u)),
        varP(utl::zeros(3u, 3u)),
        varA(utl::zeros(3u, 3u)),
        varN(utl::zeros(3u, 3u)),
        varD(utl::zeros(3u)),
        varI(utl::zeros(3u)),
        varT(utl::zeros(3u)),
        Pst(utl::zeros(3u)),
        Gst(utl::zeros(3u)),
        Qst(utl::zeros(3u)),
        Cst(utl::zeros(3u)),
        Fst(utl::zeros(3u)),
        genomescan(emptyloci(arch)),
        networkscan(emptyconnexions(arch)),
        EI(0.0),
        SI(0.0),
        RI(0.0)
    {}

    void analyze(const MetaPop&, const Param&, const GenArch&);

    // Getters called in tests
    double getEI() const
    {
        return EI;
    }
    double getSI() const
    {
        return SI;
    }
    double getRI() const
    {
        return RI;
    }


    double getVarP(const size_t &t) const // used in test
    {
        return varP[t][2u];
    }

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

    std::vector<Locus> emptyloci(const GenArch&) const;
    std::vector<Connexion> emptyconnexions(const GenArch&) const;

    std::vector<std::vector<size_t> > counts; // per habitat per ecotype

    std::vector<std::vector<std::vector<double> > > means; // per trait per habitat per ecotype

    std::vector<std::vector<double> > varG; // per trait per ecotype
    std::vector<std::vector<double> > varP; // per trait per ecotype
    std::vector<std::vector<double> > varA; // per trait per ecotype
    std::vector<std::vector<double> > varN; // per trait per ecotype
    std::vector<double> varD; // per trait
    std::vector<double> varI; // per trait
    std::vector<double> varT; // per trait

    std::vector<double> Pst; // per trait
    std::vector<double> Gst; // per trait
    std::vector<double> Qst; // per trait
    std::vector<double> Cst; // per trait
    std::vector<double> Fst; // per trait

    std::vector<Locus> genomescan; // per locus
    std::vector<Connexion> networkscan; // per edge

    double EI;
    double SI;
    double RI;

};

double Xst(const std::vector<double>&, const std::vector<size_t>&);

#endif
