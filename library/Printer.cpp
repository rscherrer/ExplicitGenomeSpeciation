#include "Printer.h"

// Constructor
Printer::Printer(const std::string &orderfile, const bool &datsave) :
    filenames(whattosave(orderfile)),
    files({ })
{

    if (datsave) {

        files.reserve(filenames.size());

        // Open files
        for (size_t f = 0u; f < filenames.size(); ++f) {

            const std::string filename = filenames[f] + ".dat";
            Stream out(new std::ofstream);
            out->open(filename.c_str(), std::ios::binary);
            if (!out->is_open()) {
                std::string msg = "Unable to open output file " + filename;
                throw std::runtime_error(msg);
            }
            files.push_back(out);
        }

    }

}

Printer::~Printer() {
    shutdown();
}

void Printer::shutdown()
{
    // Close files
    for (size_t f = 0u; f < files.size(); ++f)
        files[f]->close();
}

void Printer::print(const size_t &t, const Collector &c, const MetaPop &m)
{

    for (size_t f = 0u; f < filenames.size(); ++f) {

        if (filenames[f] == "time")
            stf::write(utl::size2dbl(t), files[f]);
        else if (filenames[f] == "population_size")
            stf::write(utl::size2dbl(c.counts[2u][2u]), files[f]);
        else if (filenames[f] == "ecotype_size")
            for (size_t eco = 0u; eco < 2u; ++eco)
                stf::write(utl::size2dbl(c.counts[2u][eco]), files[f]);
        else if (filenames[f] == "resources")
            for (size_t hab = 0u; hab < 2u; ++hab)
                for (size_t res = 0u; res < 2u; ++res)
                    stf::write(m.resources[hab][res], files[f]);
        else if (filenames[f] == "means")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.means[trait][2u][2u], files[f]);
        else if (filenames[f] == "ecotype_means")
            for (size_t trait = 0u; trait < 3u; ++trait)
                for (size_t eco = 0u; eco < 2u; ++eco)
                    stf::write(c.means[trait][2u][eco], files[f]);
        else if (filenames[f] == "varP")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varP[trait][2u], files[f]);
        else if (filenames[f] == "varG")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varG[trait][2u], files[f]);
        else if (filenames[f] == "varA")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varA[trait][2u], files[f]);
        else if (filenames[f] == "varD")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varD[trait], files[f]);
        else if (filenames[f] == "varI")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varI[trait], files[f]);
        else if (filenames[f] == "varN")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varN[trait][2u], files[f]);
        else if (filenames[f] == "varT")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.varT[trait], files[f]);
        else if (filenames[f] == "Pst")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.Pst[trait], files[f]);
        else if (filenames[f] == "Gst")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.Gst[trait], files[f]);
        else if (filenames[f] == "Qst")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.Qst[trait], files[f]);
        else if (filenames[f] == "Cst")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.Cst[trait], files[f]);
        else if (filenames[f] == "Fst")
            for (size_t trait = 0u; trait < 3u; ++trait)
                stf::write(c.Fst[trait], files[f]);
        else if (filenames[f] == "EI")
             stf::write(c.EI, files[f]);
        else if (filenames[f] == "SI")
             stf::write(c.SI, files[f]);
        else if (filenames[f] == "RI")
             stf::write(c.RI, files[f]);
        else if (filenames[f] == "genome_varP")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varP[2u], files[f]);
        else if (filenames[f] == "genome_varG")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varG[2u], files[f]);
        else if (filenames[f] == "genome_varA")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varA[2u], files[f]);
        else if (filenames[f] == "genome_varD")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varD, files[f]);
        else if (filenames[f] == "genome_varI")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varI, files[f]);
        else if (filenames[f] == "genome_varN")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].varN[2u], files[f]);
        else if (filenames[f] == "genome_Pst")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].Pst, files[f]);
        else if (filenames[f] == "genome_Gst")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].Gst, files[f]);
        else if (filenames[f] == "genome_Qst")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].Qst, files[f]);
        else if (filenames[f] == "genome_Cst")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].Cst, files[f]);
        else if (filenames[f] == "genome_Fst")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].Fst, files[f]);
        else if (filenames[f] == "genome_alpha")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].alpha, files[f]);
        else if (filenames[f] == "genome_meang")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].meanG, files[f]);
        else if (filenames[f] == "genome_freq")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                stf::write(c.genomescan[l].freqs[2u], files[f]);
        else if (filenames[f] == "genome_freqs")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                for (size_t e = 0u; e < 2u; ++e)
                    stf::write(c.genomescan[l].freqs[e], files[f]);
        else if (filenames[f] == "genome_hobs")
            for (size_t l = 0u; l < c.genomescan.size(); ++l)
                for (size_t e = 0u; e < 2u; ++e)
                    stf::write(c.genomescan[l].hobs[e], files[f]);
        else if (filenames[f] == "network_corgen")
            for (size_t e = 0u; e < c.networkscan.size(); ++e)
                stf::write(c.networkscan[e].corgen, files[f]);
        else if (filenames[f] == "network_corbreed")
            for (size_t e = 0u; e < c.networkscan.size(); ++e)
                stf::write(c.networkscan[e].corbreed, files[f]);
        else if (filenames[f] == "network_corfreq")
            for (size_t e = 0u; e < c.networkscan.size(); ++e)
                stf::write(c.networkscan[e].corfreq, files[f]);
        else if (filenames[f] == "network_avgi")
            for (size_t e = 0u; e < c.networkscan.size(); ++e)
                stf::write(c.networkscan[e].avgi, files[f]);
        else if (filenames[f] == "network_avgj")
            for (size_t e = 0u; e < c.networkscan.size(); ++e)
                stf::write(c.networkscan[e].avgj, files[f]);
        else if (filenames[f] == "individual_ecotype")
            for (size_t i = 0u; i < m.getSize(); ++i)
                stf::write(utl::size2dbl(m.getEcotype(i)), files[f]);
        else if (filenames[f] == "individual_habitat")
            for (size_t i = 0u; i < m.getSize(); ++i)
                stf::write(utl::size2dbl(m.getHabitat(i)), files[f]);
        else if (filenames[f] == "individual_trait")
            for (size_t i = 0u; i < m.getSize(); ++i)
                for (size_t trait = 0u; trait < 3u; ++trait)
                    stf::write(m.getTrait(i, trait), files[f]);
        else if (filenames[f] == "individual_midparent")
            for (size_t i = 0u; i < m.getSize(); ++i)
                for (size_t trait = 0u; trait < 3u; ++trait)
                    stf::write(m.getMidparent(i, trait), files[f]);
    }
}

std::vector<std::string> Printer::whattosave(const std::string &filename) const
{
    if (filename == "") {

        // Save the following variables if none defined

        return {

            "time",
            "population_size",
            "ecotype_size", // per ecotype
            "resources", // per habitat per resource
            "means", // per trait per ecotype
            "ecotype_means", // per trait per ecotype
            "varP", // per trait
            "varG", // per trait
            "varA", // per trait
            "varD", // per trait
            "varI", // per trait
            "varN", // per trait
            "varT", // per trait
            "Pst", // per trait
            "Gst", // per trait
            "Qst", // per trait
            "Cst", // per trait
            "Fst", // per trait
            "EI",
            "SI",
            "RI",
            "genome_varP", // per locus
            "genome_varG", // per locus
            "genome_varA", // per locus
            "genome_varD", // per locus
            "genome_varI", // per locus
            "genome_varN", // per locus
            "genome_Pst", // per locus
            "genome_Gst", // per locus
            "genome_Qst", // per locus
            "genome_Cst", // per locus
            "genome_Fst", // per locus
            "genome_alpha", // per locus
            "genome_meang", // per locus
            "genome_freq", // per locus
            "genome_freqs", // per locus per ecotype
            "genome_hobs", // per locus per ecotype
            "network_corgen", // per edge
            "network_corbreed", // per edge
            "network_corfreq", // per edge
            "network_avgi", // per edge
            "network_avgj", // per edge
            "individual_ecotype", // per individual
            "individual_habitat", // per individual
            "individual_trait", // per individual per trait
            "individual_midparent" // per individual per trait

        };

    }

    std::vector<std::string> variables;

    // Or read file defining what variables to save

    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        std::string msg = "Unable to open save file ";
        throw std::runtime_error(msg + filename);
    }

    std::string input;
    while (file >> input) variables.push_back(input);
    file.close();
    return variables;

}

