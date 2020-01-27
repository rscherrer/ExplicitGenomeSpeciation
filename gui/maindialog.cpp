#include "gui/maindialog.h"
#include "ui_maindialog.h"
#include "library/Simul.h"
#include <sstream>
#include <QThread>
#include <QColor>

MainDialog::MainDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::MainDialog)
{

  ui->setupUi(this);
  setup_spinboxes();

  ui->plot->addGraph();
  ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssDisc);
  ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone);
  ui->plot->graph(0)->setPen(QPen(Qt::red));

  ui->plot->xAxis->grid()->setVisible(false);
  ui->plot->yAxis->grid()->setVisible(false);

  QCPPlotTitle *fst_title = new QCPPlotTitle(ui->plot, "Fst");
  ui->plot->plotLayout()->insertRow(0);
  ui->plot->plotLayout()->addElement(0, 0, fst_title);
  ui->plot->xAxis->setLabel("Locus");
  ui->plot->yAxis->setLabel("Fst");

  ui->plot->xAxis->setRange(0,90);
  ui->plot->yAxis->setRange(0,1);


  // overall population size graph
  ui->plot_popsize->addGraph();
  ui->plot_popsize->graph(0)->setScatterStyle(QCPScatterStyle::ssDisc);
  ui->plot_popsize->graph(0)->setPen(QPen(Qt::red));
  ui->plot_popsize->graph(0)->setName("Deme 0");

  // deme 1 population size graph
  ui->plot_popsize->addGraph(ui->plot_popsize->xAxis, ui->plot_popsize->yAxis);
  ui->plot_popsize->graph(1)->setScatterStyle(QCPScatterStyle::ssDisc);
  ui->plot_popsize->graph(1)->setPen(QPen(Qt::blue));
  ui->plot_popsize->graph(1)->setName("Deme 1");

  ui->plot_popsize->yAxis2->setLabel("Deme Size");
  // ui->plot_popsize->yAxis2->setVisible(true);


  QCPPlotTitle *popsize_title = new QCPPlotTitle(ui->plot_popsize, "Population Size");
  ui->plot_popsize->plotLayout()->insertRow(0);
  ui->plot_popsize->plotLayout()->addElement(0, 0, popsize_title);
  ui->plot_popsize->xAxis->setLabel("Time");
  ui->plot_popsize->yAxis->setLabel("Size Deme 0");
  ui->plot_popsize->yAxis2->setLabel("Size Deme 1");

   


  ui->plot_eco_trait->addGraph();
  ecoBars_0 = new QCPBars(ui->plot_eco_trait->xAxis,
                         ui->plot_eco_trait->yAxis);
  ecoBars_0->setName("Eco Trait deme 0");
  ecoBars_0->setPen(QPen(Qt::red));
  ecoBars_0->setBrush(QBrush(QColor(255,0,0, static_cast<int>(0.5 * 255))));
  ecoBars_0->setAntialiased(false);
  ecoBars_0->setAntialiasedFill(false);

  ui->plot_eco_trait->addGraph();
  ecoBars_1 = new QCPBars(ui->plot_eco_trait->xAxis,
                         ui->plot_eco_trait->yAxis);
  ecoBars_1->setName("Eco Trait deme 1");
  ecoBars_1->setPen(QPen(Qt::blue));
  ecoBars_1->setBrush(QBrush(QColor(0,0,255, static_cast<int>(0.5 * 255))));
  ecoBars_1->setAntialiased(false);
  ecoBars_1->setAntialiasedFill(false);


  QCPPlotTitle *eco_title = new QCPPlotTitle(ui->plot_eco_trait, "Ecological Trait");
  ui->plot_eco_trait->plotLayout()->insertRow(0);
  ui->plot_eco_trait->plotLayout()->addElement(0, 0, eco_title);
  ui->plot_eco_trait->xAxis->setLabel("Trait Value");
  ui->plot_eco_trait->yAxis->setLabel("Count");


  ui->plot_sex_trait->addGraph();
  sexBars_0 = new QCPBars(ui->plot_sex_trait->xAxis,
                          ui->plot_sex_trait->yAxis);
  sexBars_0->setName("Sex Trait Deme 0");
  sexBars_0->setPen(QPen(Qt::red));
  sexBars_0->setBrush(QColor(255,0,0, static_cast<int>(0.5 * 255)));
  sexBars_0->setAntialiased(false);
  sexBars_0->setAntialiasedFill(false);


  ui->plot_sex_trait->addGraph();
  sexBars_1 = new QCPBars(ui->plot_sex_trait->xAxis,
                          ui->plot_sex_trait->yAxis);
  sexBars_1->setName("Sex Trait Deme 1");
  sexBars_1->setPen(QPen(Qt::blue));
  sexBars_1->setBrush(QColor(0,0,255, static_cast<int>(0.5 * 255)));
  sexBars_1->setAntialiased(false);
  sexBars_1->setAntialiasedFill(false);



  QCPPlotTitle *sex_title = new QCPPlotTitle(ui->plot_sex_trait, "Mating Preference");
  ui->plot_sex_trait->plotLayout()->insertRow(0);
  ui->plot_sex_trait->plotLayout()->addElement(0, 0, sex_title);
  ui->plot_sex_trait->xAxis->setLabel("Mating Preference");
  ui->plot_sex_trait->yAxis->setLabel("Count");

  ui->plot_neu_trait->addGraph();

  neuBars_0 = new QCPBars(ui->plot_neu_trait->xAxis,
                        ui->plot_neu_trait->yAxis);
  neuBars_0->setName("Neutral Trait Deme 0");
  neuBars_0->setPen(QPen(Qt::red));
  neuBars_0->setBrush(QColor(255, 0, 0, static_cast<int>(0.5 * 255)));
  neuBars_0->setAntialiased(false);
  neuBars_0->setAntialiasedFill(false);

  neuBars_1 = new QCPBars(ui->plot_neu_trait->xAxis,
                        ui->plot_neu_trait->yAxis);
  neuBars_1->setName("Neutral Trait Deme 0");
  neuBars_1->setPen(QPen(Qt::blue));
  neuBars_1->setBrush(QColor(0, 0, 255, static_cast<int>(0.5 * 255)));
  neuBars_1->setAntialiased(false);
  neuBars_1->setAntialiasedFill(false);

  QCPPlotTitle *neu_title = new QCPPlotTitle(ui->plot_neu_trait, "Neutral Trait");
  ui->plot_neu_trait->plotLayout()->insertRow(0);
  ui->plot_neu_trait->plotLayout()->addElement(0, 0, neu_title);
  ui->plot_neu_trait->xAxis->setLabel("Trait Value");
  ui->plot_neu_trait->yAxis->setLabel("Count");
}

MainDialog::~MainDialog()
{
    delete ui;
}

void MainDialog::update_plot_popsize(int t,
                                     size_t n_0,
                                     size_t n_1) {
    pop_x.append(t);

    pop_y_0.append(n_0);
    pop_y_1.append(n_1);

    ui->plot_popsize->graph(0)->clearData();
    ui->plot_popsize->graph(0)->setData(pop_x, pop_y_0);

    ui->plot_popsize->graph(1)->clearData();
    ui->plot_popsize->graph(1)->setData(pop_x, pop_y_1);

    ui->plot_popsize->rescaleAxes();
    ui->plot_popsize->replot();
    ui->plot_popsize->update();

}

void MainDialog::plot_fst(const std::vector<double>& v)
{

    QVector<double> fst_x(static_cast<int>(v.size()));
    std::iota(fst_x.begin(), fst_x.end(), 1);

    QVector<double> fst_y = QVector<double>::fromStdVector(v);

    ui->plot->graph(0)->clearData();
    ui->plot->graph(0)->setData(fst_x, fst_y);
    ui->plot->rescaleAxes();
    ui->plot->replot();
    ui->plot->update();
}


void update_barplot(const std::vector< double>& v,
                    QCPBars * barplot,
                    double& max_y_val,
                    double& min_x,
                    double& max_x) {

    auto min_max = std::minmax_element(v.begin(), v.end());

    if(*min_max.first < min_x) min_x = *min_max.first;
    if(*min_max.second > max_x) max_x = *min_max.second;

    int num_bins = 30;
    double bin_size = 1.0 * (*min_max.second - *min_max.first) / num_bins;

    barplot->setWidth(bin_size);

    std::vector<int> bins(30, 0);

    for(size_t i = 0; i < 30; ++i) {
        double left = *min_max.first + i * bin_size;
        double right = *min_max.first + (i+1) * bin_size;
        for(auto it = v.begin(); it != v.end(); ++it) {
            if( (*it) >= left && (*it) < right) {
                bins[ i ]++;
            }
        }
    }

    QVector<double> xvals;
    QVector<double> yvals;

    for(size_t i = 0; i < bins.size(); ++i) {
        xvals.append(*min_max.first + i * bin_size);
        yvals.append(bins[i]);
        if(bins[i] > max_y_val) max_y_val = bins[i];
    }

    barplot->clearData();
    barplot->setData(xvals, yvals);
}


void plot_barplot(QCustomPlot* UI,
                  QCPBars * barplot_0,
                  QCPBars * barplot_1,
                  const std::vector<double>& v_0,
                  const std::vector<double>& v_1) {

    double max_y_val = -1;
    double min_x = 1e6;
    double max_x = -1e6;

    update_barplot(v_0, barplot_0, max_y_val, min_x, max_x);
    update_barplot(v_1, barplot_1, max_y_val, min_x, max_x);


    UI->xAxis->setRange(0.9 * min_x, 1.1 * max_x);

    UI->yAxis->setRange(0, max_y_val * 1.05);
    UI->yAxis->setRange(0, max_y_val * 1.05);

    UI->replot();
    UI->update();

    return;
}

void MainDialog::on_run_button_clicked()
{
    is_running = true;
    try
    {
        // clean up old data
        pop_x.clear();
        pop_y_0.clear();
        pop_y_1.clear();

        // Create the parameters from the GUI
        Param pars = createPars();

        // Random number generator
        rnd::rng.seed(pars.seed);

        // Create a genetic architecture
        GenArch arch = GenArch(pars);

        // Create a metapopulation with two demes
        MetaPop metapop = MetaPop(pars, arch);

        // Create an analytical module
        Collector collector = Collector(arch);

        // Loop through time
        for (int t = -pars.tburnin; t < pars.tend; ++t) {
            if(!is_running) {
                break;
            }

            // flag to make sure the GUI is updated realtime
            QApplication::processEvents();

            if (t == 0) metapop.exitburnin();

            // Life cycle of the metapopulation
            metapop.cycle(pars, arch);

            // Is the population still there?
            if (metapop.isextinct()) {
                std::cout << "The population went extinct at t = " << t << '\n';
                break;
            }

            // Analyze the metapopulation if needed
            if (timetosave(t, pars)) {

                std::stringstream s;
                s << "t: " << t;
                ui->output->appendPlainText(QString::fromStdString(s.str()));

                collector.analyze(metapop, pars, arch);

                std::vector<double> fst_vals = collector.get_Fst();
                plot_fst(fst_vals);

                std::vector<double> eco_trait_0 = collector.get_eco_trait_deme(metapop, 0);
                std::vector<double> eco_trait_1 = collector.get_eco_trait_deme(metapop, 1);

                std::vector<double> sex_trait_0 = collector.get_sex_trait_deme(metapop, 0);
                std::vector<double> sex_trait_1 = collector.get_sex_trait_deme(metapop, 1);

                std::vector<double> neu_trait_0 = collector.get_neu_trait_deme(metapop, 0);
                std::vector<double> neu_trait_1 = collector.get_neu_trait_deme(metapop, 1);

                //plot_hist(ui->plot_eco_trait, ecoBars, eco_trait_vals);
                plot_barplot(ui->plot_eco_trait, ecoBars_0, ecoBars_1,
                             eco_trait_0, eco_trait_1);

                plot_barplot(ui->plot_sex_trait, sexBars_0, sexBars_1,
                             sex_trait_0, sex_trait_1);

                plot_barplot(ui->plot_neu_trait, neuBars_0, neuBars_1,
                             neu_trait_0, neu_trait_1);


                update_plot_popsize(t,
                                    metapop.getDemeSize(0),
                                    metapop.getDemeSize(1));
            }
        }
        // Show output
        {
            std::stringstream s;
            s << "Done\n";
            ui->output->appendPlainText(QString::fromStdString(s.str()));
        }
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
    }
}

void MainDialog::on_pushButton_clicked()
{
    // boolean flag that stops simulation
    is_running = false;
}



void MainDialog::setup_spinboxes() {
    // set spinboxes to default values,
    // where default values come from Params.h:

    Param temp_pars;
    ui->box_rdynamics->setValue(static_cast<int>(temp_pars.rdynamics));
    ui->box_trenewal->setValue(temp_pars.trenewal);
    ui->box_capacity->setValue(temp_pars.capacity);
    ui->box_replenish->setValue(temp_pars.replenish);
    ui->box_hsymmetry->setValue(temp_pars.hsymmetry);
    ui->box_ecosel->setValue(temp_pars.ecosel);
    ui->box_dispersal->setValue(temp_pars.dispersal);
    ui->box_birth->setValue(temp_pars.birth);
    ui->box_survival->setValue(temp_pars.survival);
    ui->box_sexsel->setValue(temp_pars.sexsel);
    ui->box_matingcost->setValue(temp_pars.matingcost);
    ui->box_maxfeed->setValue(temp_pars.maxfeed);
    ui->box_nloci->setValue(static_cast<int>(temp_pars.nloci));
    ui->box_nchrom->setValue(static_cast<int>(temp_pars.nchrom));
    ui->box_mutation->setValue(temp_pars.mutation);
    ui->box_recombination->setValue(temp_pars.recombination);
    ui->box_allfreq->setValue(temp_pars.allfreq);
    ui->box_effectshape->setValue(temp_pars.effectshape);
    ui->box_effectscale->setValue(temp_pars.effectscale);
    ui->box_interactionshape->setValue(temp_pars.interactionshape);
    ui->box_interactionscale->setValue(temp_pars.interactionscale);
    ui->box_dominancevar->setValue(temp_pars.dominancevar);
    ui->box_tburnin->setValue(temp_pars.tburnin);
    ui->box_tend->setValue(temp_pars.tend);
    ui->box_tsave->setValue(temp_pars.tsave);
    ui->box_record->setValue(temp_pars.record);
    ui->box_seed->setValue(static_cast<int>(temp_pars.seed));
    ui->box_ntrials->setValue(static_cast<int>(temp_pars.ntrials));


    ui->box_demesizes_0->setValue(static_cast<int>(temp_pars.demesizes[0]));
    ui->box_demesizes_1->setValue(static_cast<int>(temp_pars.demesizes[1]));

    ui->box_nvertices_0->setValue(static_cast<int>(temp_pars.nvertices[0]));
    ui->box_nvertices_1->setValue(static_cast<int>(temp_pars.nvertices[1]));
    ui->box_nvertices_2->setValue(static_cast<int>(temp_pars.nvertices[2]));

    ui->box_nedges_0->setValue(static_cast<int>(temp_pars.nedges[0]));
    ui->box_nedges_1->setValue(static_cast<int>(temp_pars.nedges[1]));
    ui->box_nedges_2->setValue(static_cast<int>(temp_pars.nedges[2]));

    ui->box_scaleA_0->setValue(temp_pars.scaleA[0]);
    ui->box_scaleA_1->setValue(temp_pars.scaleA[1]);
    ui->box_scaleA_2->setValue(temp_pars.scaleA[2]);

    ui->box_scaleD_0->setValue(temp_pars.scaleD[0]);
    ui->box_scaleD_1->setValue(temp_pars.scaleD[1]);
    ui->box_scaleD_2->setValue(temp_pars.scaleD[2]);

    ui->box_scaleI_0->setValue(temp_pars.scaleI[0]);
    ui->box_scaleI_1->setValue(temp_pars.scaleI[1]);
    ui->box_scaleI_2->setValue(temp_pars.scaleI[2]);

    ui->box_scaleE_0->setValue(temp_pars.scaleE[0]);
    ui->box_scaleE_1->setValue(temp_pars.scaleE[1]);
    ui->box_scaleE_2->setValue(temp_pars.scaleE[2]);

    ui->box_locusE_0->setValue(temp_pars.locusE[0]);
    ui->box_locusE_1->setValue(temp_pars.locusE[1]);
    ui->box_locusE_2->setValue(temp_pars.locusE[2]);

    ui->box_skews_0->setValue(temp_pars.skews[0]);
    ui->box_skews_0->setValue(temp_pars.skews[1]);
    ui->box_skews_0->setValue(temp_pars.skews[2]);
}

Param MainDialog::createPars()
{
    // read parameters from GUI
    Param pars;
    pars.rdynamics          = static_cast<size_t>(ui->box_rdynamics->value());
    pars.trenewal           = ui->box_trenewal->value();
    pars.capacity           = ui->box_capacity->value();
    pars.replenish          = ui->box_replenish->value();
    pars.hsymmetry          = ui->box_hsymmetry->value();
    pars.ecosel             = ui->box_ecosel->value();
    pars.dispersal          = ui->box_dispersal->value();
    pars.birth              = ui->box_birth->value();
    pars.survival           = ui->box_survival->value();
    pars.sexsel             = ui->box_sexsel->value();
    pars.matingcost         = ui->box_matingcost->value();
    pars.maxfeed            = ui->box_maxfeed->value();
    pars.nloci              = static_cast<size_t>(ui->box_nloci->value());
    pars.nchrom             = static_cast<size_t>(ui->box_nchrom->value());
    pars.mutation           = ui->box_mutation->value();
    pars.recombination      = ui->box_recombination->value();
    pars.allfreq            = ui->box_allfreq->value();
    pars.effectshape        = ui->box_effectshape->value();
    pars.effectscale        = ui->box_effectscale->value();
    pars.interactionshape   = ui->box_interactionshape->value();
    pars.interactionscale   = ui->box_interactionscale->value();
    pars.dominancevar       = ui->box_dominancevar->value();
    pars.tburnin            = static_cast<int>(ui->box_tburnin->value());
    pars.tend               = static_cast<int>(ui->box_tend->value());
    pars.tsave              = static_cast<int>(ui->box_tsave->value());
    pars.record             = static_cast<bool>(ui->box_record->value());
    pars.seed               = static_cast<size_t>(ui->box_seed->value());
    pars.ntrials            = static_cast<size_t>(ui->box_ntrials->value());

    pars.demesizes[0]       = static_cast<size_t>(ui->box_demesizes_0->value());
    pars.demesizes[1]       = static_cast<size_t>(ui->box_demesizes_1->value());

    pars.nvertices[0]       = static_cast<size_t>(ui->box_nvertices_0->value());
    pars.nvertices[1]       = static_cast<size_t>(ui->box_nvertices_1->value());
    pars.nvertices[2]       = static_cast<size_t>(ui->box_nvertices_2->value());

    pars.nedges[0]          = static_cast<size_t>(ui->box_nedges_0->value());
    pars.nedges[1]          = static_cast<size_t>(ui->box_nedges_1->value());
    pars.nedges[2]          = static_cast<size_t>(ui->box_nedges_2->value());

    pars.scaleA[0]          = ui->box_scaleA_0->value();
    pars.scaleA[1]          = ui->box_scaleA_1->value();
    pars.scaleA[2]          = ui->box_scaleA_2->value();

    pars.scaleD[0]          = ui->box_scaleD_0->value();
    pars.scaleD[1]          = ui->box_scaleD_1->value();
    pars.scaleD[2]          = ui->box_scaleD_2->value();

    pars.scaleI[0]          = ui->box_scaleI_0->value();
    pars.scaleI[1]          = ui->box_scaleI_1->value();
    pars.scaleI[2]          = ui->box_scaleI_2->value();

    pars.scaleE[0]          = ui->box_scaleE_0->value();
    pars.scaleE[1]          = ui->box_scaleE_1->value();
    pars.scaleE[2]          = ui->box_scaleE_2->value();

    pars.locusE[0]          = ui->box_locusE_0->value();
    pars.locusE[1]          = ui->box_locusE_1->value();
    pars.locusE[2]          = ui->box_locusE_2->value();

    pars.skews[0]           = ui->box_skews_0->value();
    pars.skews[1]           = ui->box_skews_1->value();
    pars.skews[2]           = ui->box_skews_2->value();



    pars.seed = static_cast<size_t>(ui->rng_seed->value());
    pars.tend = static_cast<int>(ui->num_gen->value());

    std::stringstream s_p;
    s_p << "parameters read from GUI\n";
    s_p << "With the following values:\n";

    // display parameters:

    s_p << "rdynamics: "            << pars.rdynamics << "\n";
    s_p << "trenewal: "             << pars.trenewal << "\n";
    s_p << "capacity: "             << pars.capacity << "\n";
    s_p << "replenish: "            << pars.replenish << "\n";
    s_p << "hsymmetry: "            << pars.hsymmetry << "\n";
    s_p << "ecosel: "               << pars.ecosel << "\n";
    s_p << "dispersal: "            << pars.dispersal << "\n";
    s_p << "birth: "                << pars.birth << "\n";
    s_p << "survival: "             << pars.survival << "\n";
    s_p << "sexsel: "               << pars.sexsel << "\n";
    s_p << "matingcost: "           << pars.matingcost << "\n";
    s_p << "maxfeed: "              << pars.maxfeed << "\n";
    s_p << "demesizes: "            << "{" << pars.demesizes[0] << "," << pars.demesizes[1] << "}\n";
    s_p << "nloci: "                << pars.nloci << "\n";
    s_p << "nvertices: "            << "{" << pars.nvertices[0] << "," << pars.nvertices[1] << "," << pars.nvertices[2] << "}\n";
    s_p << "nedges: "               << "{" << pars.nedges[0]    << "," << pars.nedges[1]    << "," << pars.nedges[2] << "}\n";
    s_p << "nchrom: "               << pars.nchrom << "\n";
    s_p << "mutation: "             << pars.mutation << "\n";
    s_p << "recombination: "        << pars.recombination << "\n";
    s_p << "allfreq: "              << pars.allfreq << "\n";
    s_p << "scaleA: "               << "{" << pars.scaleA[0] << "," << pars.scaleA[1] << "," << pars.scaleA[2] << "}\n";
    s_p << "scaleD: "               << "{" << pars.scaleD[0] << "," << pars.scaleD[1] << "," << pars.scaleD[2] << "}\n";
    s_p << "scaleI: "               << "{" << pars.scaleI[0] << "," << pars.scaleI[1] << "," << pars.scaleI[2] << "}\n";
    s_p << "scaleE: "               << "{" << pars.scaleE[0] << "," << pars.scaleE[1] << "," << pars.scaleE[2] << "}\n";
    s_p << "locusE: "               << "{" << pars.locusE[0] << "," << pars.locusE[1] << "," << pars.locusE[2] << "}\n";
    s_p << "skews: "                << "{" << pars.skews[0]  << "," << pars.skews[1]  << "," << pars.skews[2]  << "}\n";
    s_p << "effectshape: "          << pars.effectshape << "\n";
    s_p << "effectscale: "          << pars.effectscale << "\n";
    s_p << "interactionshape: "     << pars.interactionshape << "\n";
    s_p << "interactionscale: "     << pars.interactionscale << "\n";
    s_p << "dominancevar: "         << pars.dominancevar << "\n";
    s_p << "tburnin: "              << pars.tburnin << "\n";
    s_p << "tend: "                 << pars.tend << "\n";
    s_p << "tsave: "                << pars.tsave << "\n";
    s_p << "record: "               << pars.record << "\n";
    s_p << "seed: "                 << pars.seed << "\n";
    s_p << "ntrials: "              << pars.ntrials << "\n";

    ui->output->appendPlainText(QString::fromStdString(s_p.str()));

    return pars;
}




/// unused code /////
/// (for now)   /////

void plot_hist(QCustomPlot* UI, QCPBars * barplot,
               const std::vector<double>& v) {

    auto min_max = std::minmax_element(v.begin(), v.end());

    int num_bins = 30;
    double bin_size = 1.0 * (*min_max.second - *min_max.first) / num_bins;

    barplot->setWidth(bin_size);

    std::vector<int> bins(30, 0);

    for(size_t i = 0; i < 30; ++i) {
        double left = *min_max.first + i * bin_size;
        double right = *min_max.first + (i+1) * bin_size;
        for(auto it = v.begin(); it != v.end(); ++it) {
            if( (*it) >= left && (*it) < right) {
                bins[ i ]++;
            }
        }
    }

    QVector<double> xvals;
    QVector<double> yvals;

    double max_y_val = -1;
    for(size_t i = 0; i < bins.size(); ++i) {
        xvals.append(*min_max.first + i * bin_size);
        yvals.append(bins[i]);
        if(bins[i] > max_y_val) max_y_val = bins[i];
    }

    barplot->clearData();

    UI->xAxis->setRange(0.9 * *min_max.first,
                                        1.1 * *min_max.second);

    UI->yAxis->setRange(0, max_y_val * 1.05);

    barplot->setData(xvals, yvals);
    UI->replot();
    UI->update();

    return;
}
