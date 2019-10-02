#ifndef EXPLICITGENOMESPECIATION_STATS_H
#define EXPLICITGENOMESPECIATION_STATS_H

#include "Types.h"
#include "Utilities.h"
#include "Deme.h"
#include "Output.h"

typedef std::vector<Deme> vecPop;

class Stats
{

public:

    Stats(const GenArch &arch) :
        time(0u),
        genomescan(makeEmptyGenomeScan(arch)),
        traitstats(makeEmptyTraits()),
        resources(utl::matzeros(2u, 2u)),
        popcounts(utl::uzeros(2u)),
        nfemales(utl::uzeros(2u)),
        ecocounts(utl::uzeros(2u)),
        totcount(0u),
        EI(0.0),
        SI(0.0),
        RI(0.0)
    {}

    class Locus
    {

        friend class Stats;

    public:

        Locus(const size_t &l, const size_t &envar) :
            id(l),
            genocounts(utl::uzeros(3u)),
            genosumgens(utl::zeros(3u)),
            genossqgens(utl::zeros(3u)),
            genobetas(utl::zeros(3u)),
            genoexpecs(utl::zeros(3u)),
            genodeltas(utl::zeros(3u)),
            ecostats(makeEmptyEcotypes()),
            sumgen(0.0),
            ssqgen(0.0),
            p(0.0),
            meanq(0.0),
            varq(0.0),
            covqg(0.0),
            alpha(0.0),
            varE(envar),
            varP(0.0),
            varG(0.0),
            varA(0.0),
            varD(0.0),
            varI(0.0),
            varN(0.0),
            varS(0.0),
            varT(0.0),
            Pst(0.0),
            Gst(0.0),
            Qst(0.0),
            Cst(0.0),
            Fst(0.0)
        {}

        class Ecotype
        {

            friend class Stats;
            friend class Stats::Locus;

        public:

            Ecotype(const size_t &e) :
                eco(e),
                sumgen(0.0),
                ssqgen(0.0),
                p(0.0),
                genocounts(utl::uzeros(3u)),
                genosumgens(utl::zeros(3u)),
                genossqgens(utl::zeros(3u)),
                varP(0.0),
                varG(0.0),
                varA(0.0),
                varN(0.0)
            {}

            void setup();
            void setVarG();
            void setVarP(const double&);
            void setVarA(const vecDbl&);
            void setVarN(const vecUns&, const vecDbl&);

        private:

            size_t eco;

            double sumgen;
            double ssqgen;
            double p;

            vecUns genocounts;
            vecDbl genosumgens;
            vecDbl genossqgens;

            double varP;
            double varG;
            double varA;
            double varN;

        };

        double getVarP(const size_t &e) const { return ecostats[e].varP; }
        double getVarG(const size_t &e) const { return ecostats[e].varG; }
        double getVarA(const size_t &e) const { return ecostats[e].varA; }
        double getVarN(const size_t &e) const { return ecostats[e].varN; }
        double getP(const size_t &e) const { return ecostats[e].p; }

        void accumulate(const vecPop&);
        void regress();
        void setEcotypes();
        void setVarG();
        void setVarP();
        void setVarA();
        void setVarD();
        void setVarI();
        void setVarN();
        void setVarS();
        void setVarT();
        void setPst();
        void setGst();
        void setQst();
        void setCst();
        void setFst();

        void setTotCount(const double&);
        void setEcoCounts(const vecUns&);

    private:

        size_t id;

        std::vector<Stats::Locus::Ecotype> makeEmptyEcotypes();

        vecUns genocounts;
        vecDbl genosumgens;
        vecDbl genossqgens;

        vecDbl genobetas;
        vecDbl genoexpecs;
        vecDbl genodeltas;

        std::vector<Stats::Locus::Ecotype> ecostats;

        double sumgen;
        double ssqgen;

        double p;
        double meanq;
        double varq;
        double covqg;
        double alpha;

        double varE;
        double varP;
        double varG;
        double varA;
        double varD;
        double varI;
        double varN;
        double varS;
        double varT;
        double Pst;
        double Gst;
        double Qst;
        double Cst;
        double Fst;

        static constexpr size_t aa = 0u;
        static constexpr size_t Aa = 1u;
        static constexpr size_t AA = 2u;

        static size_t totcount;
        static vecUns ecocounts;

    };

    class Trait
    {

        friend class Stats;

    public:

        Trait() :
            ecostats(makeEmptyEcotypes()),
            sumgen(0.0),
            sumphe(0.0),
            ssqgen(0.0),
            ssqphe(0.0),
            varP(0.0),
            varG(0.0),
            varA(0.0),
            varD(0.0),
            varI(0.0),
            varN(0.0),
            varS(0.0),
            varT(0.0),
            Pst(0.0),
            Gst(0.0),
            Qst(0.0),
            Cst(0.0),
            Fst(0.0)
        {}

        class Ecotype
        {

            friend class Stats;
            friend class Stats::Trait;

        public:

            Ecotype() :
                sumgen(0.0),
                sumphe(0.0),
                ssqgen(0.0),
                ssqphe(0.0),
                varP(0.0),
                varG(0.0),
                varA(0.0),
                varN(0.0)
            {}

            void setVarP(const size_t&);
            void setVarG(const size_t&);

        private:

            double sumgen;
            double sumphe;
            double ssqgen;
            double ssqphe;

            double varP;
            double varG;
            double varA;
            double varN;

        };

        double getMeanP() const;
        double getMeanP(const size_t&) const;
        double getSumGen(const size_t&) const;
        double getSumPhe(const size_t&) const;
        double getSsqGen(const size_t&) const;
        double getSsqPhe(const size_t&) const;
        double getVarP(const size_t&) const;
        double getVarG(const size_t&) const;
        double getVarA(const size_t&) const;
        double getVarN(const size_t&) const;

        void setVarG();
        void setVarP();
        void setVarG(const size_t&);
        void setVarP(const size_t&);
        void setPst();
        void setGst();
        void setQst();
        void setCst();
        void setFst();

        void setTotCount(const double&);
        void setEcoCounts(const vecUns&);

    private:

        std::vector<Stats::Trait::Ecotype> makeEmptyEcotypes();

        std::vector<Stats::Trait::Ecotype> ecostats;

        double sumgen;
        double sumphe;
        double ssqgen;
        double ssqphe;

        double varP;
        double varG;
        double varA;
        double varD;
        double varI;
        double varN;
        double varS;
        double varT;
        double Pst;
        double Gst;
        double Qst;
        double Cst;
        double Fst;

        static size_t totcount;
        static vecUns ecocounts;
    };

    // Used for testing
    double getEcoIsolation() const { return EI; }
    double getSpatialIsolation() const { return SI; }
    double getMatingIsolation() const { return RI; }
    double getPst(const size_t &t) const { return traitstats[t].Pst; }
    double getVarP(const size_t&, const size_t&) const;
    double getSsqPhe(const size_t&, const size_t&) const;
    double getSumPhe(const size_t&, const size_t&) const;
    double getEcoCount(const size_t &eco) const { return ecocounts[eco]; }

    void reset(const size_t&, const GenArch&);
    void analyze(const vecPop&, const GenArch&);
    void setEcoIsolation();
    void setSpatialIsolation(const vecPop&);
    void setMatingIsolation(const vecPop&, const double&, const double&);
    void save(Output&);

private:

    void accumulate(const vecPop&);
    void scanGenome(const vecPop&, const vecUns&);
    void setVariance();

    void contributeVarA(const double&, const size_t&);
    void contributeVarD(const double&, const size_t&);
    void contributeVarI(const double&, const size_t&);
    void contributeVarN(const double&, const size_t&);
    void contributeVarS(const double&, const size_t&);
    void contributeVarT(const double&, const size_t&);
    void contributeEcoVarA(const double&, const size_t&, const size_t&);
    void contributeEcoVarN(const double&, const size_t&, const size_t&);

    void write(const double&, std::ofstream *&);
    void write(const vecDbl&, std::ofstream *&);

    std::vector<Stats::Trait> makeEmptyTraits();
    std::vector<Stats::Locus> makeEmptyGenomeScan(const GenArch&);

    size_t time;

    std::vector<Stats::Locus> genomescan;
    std::vector<Stats::Trait> traitstats;

    Matrix resources;
    vecUns popcounts;
    vecUns nfemales;
    vecUns ecocounts;
    size_t totcount;

    double EI;
    double SI;
    double RI;

    static constexpr size_t x = 0u;
    static constexpr size_t y = 1u;
    static constexpr size_t z = 2u;
};


#endif
