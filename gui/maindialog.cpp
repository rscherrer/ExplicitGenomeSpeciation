#include "maindialog.h"
#include "ui_maindialog.h"
#include "library/Simul.h"
#include <sstream>
#include <QThread>

MainDialog::MainDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::MainDialog)
{
  ui->setupUi(this);

  ui->plot->addGraph();
  ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
  ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone);
  ui->plot->graph(0)->setPen(QPen(Qt::red));

  ui->plot->xAxis->grid()->setVisible(false);
  ui->plot->yAxis->grid()->setVisible(false);

  QCPPlotTitle *fst_title = new QCPPlotTitle(ui->plot, "Fst");
  ui->plot->plotLayout()->insertRow(0);
  ui->plot->plotLayout()->addElement(0, 0, fst_title);
  ui->plot->xAxis->setLabel("Locus");
  ui->plot->yAxis->setLabel("Fst");

  ui->plot_gst->addGraph();
  ui->plot_gst->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
  ui->plot_gst->graph(0)->setLineStyle(QCPGraph::lsNone);
  ui->plot_gst->graph(0)->setPen(QPen(Qt::blue));

  QCPPlotTitle *gst_title = new QCPPlotTitle(ui->plot_gst, "Gst");
  ui->plot_gst->plotLayout()->insertRow(0);
  ui->plot_gst->plotLayout()->addElement(0, 0, gst_title);
  ui->plot_gst->xAxis->setLabel("Locus");
  ui->plot_gst->yAxis->setLabel("Gst");


  ui->plot->xAxis->setRange(0,90);
  ui->plot->yAxis->setRange(0,1);
  ui->plot_gst->xAxis->setRange(0,90);
  ui->plot_gst->yAxis->setRange(0,1);


  ui->plot_eco_trait->addGraph();

  ecoBars = new QCPBars(ui->plot_eco_trait->xAxis,
                        ui->plot_eco_trait->yAxis);
  ecoBars->setName("Eco Trait");
  ecoBars->setPen(QPen(Qt::black));
  ecoBars->setAntialiased(false);
  ecoBars->setAntialiasedFill(false);

  QCPPlotTitle *eco_title = new QCPPlotTitle(ui->plot_eco_trait, "Ecological Trait");
  ui->plot_eco_trait->plotLayout()->insertRow(0);
  ui->plot_eco_trait->plotLayout()->addElement(0, 0, eco_title);
  ui->plot_eco_trait->xAxis->setLabel("Trait Value");
  ui->plot_eco_trait->yAxis->setLabel("Count");


  ui->plot_sex_trait->addGraph();
  sexBars = new QCPBars(ui->plot_sex_trait->xAxis,
                        ui->plot_sex_trait->yAxis);
  sexBars->setName("Sex Trait");
  sexBars->setPen(QPen(Qt::black));
  sexBars->setAntialiased(false);
  sexBars->setAntialiasedFill(false);

  QCPPlotTitle *sex_title = new QCPPlotTitle(ui->plot_sex_trait, "Mating Preference");
  ui->plot_sex_trait->plotLayout()->insertRow(0);
  ui->plot_sex_trait->plotLayout()->addElement(0, 0, sex_title);
  ui->plot_sex_trait->xAxis->setLabel("Mating Preference");
  ui->plot_sex_trait->yAxis->setLabel("Count");

  ui->plot_neu_trait->addGraph();

  neuBars = new QCPBars(ui->plot_neu_trait->xAxis,
                        ui->plot_neu_trait->yAxis);
  neuBars->setName("Neutral Trait");
  neuBars->setPen(QPen(Qt::black));
  neuBars->setAntialiased(false);
  neuBars->setAntialiasedFill(false);

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


void MainDialog::plot_fst(const std::vector<double>& v)
{
    fst_x.clear();
    fst_y.clear();
    for(size_t i = 0; i < v.size(); ++i) {
        fst_x.append(i);
        fst_y.append(v[i]);
    }

    ui->plot->graph(0)->clearData();
    ui->plot->graph(0)->setData(fst_x, fst_y);
    ui->plot->rescaleAxes();
    ui->plot->replot();
    ui->plot->update();
}

void MainDialog::plot_gst(const std::vector<double>& v)
{
    gst_x.clear();
    gst_y.clear();
    for(size_t i = 0; i < v.size(); ++i) {
        gst_x.append(i);
        gst_y.append(v[i]);
    }
    ui->plot_gst->graph(0)->clearData();
    ui->plot_gst->graph(0)->setData(gst_x, gst_y);
    ui->plot_gst->rescaleAxes();
    ui->plot_gst->replot();
    ui->plot_gst->update();
}

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




Param MainDialog::createPars()
{
    Param pars;
    pars.seed = static_cast<size_t>(ui->rng_seed->value());
    return pars;
}

void MainDialog::on_run_button_clicked()
{
    is_running = true;
    try
    {
        // Create the parameters from the GUI
        Param pars = createPars();

        //Show params in output
        {
            std::stringstream s;
            s << "seed: " << pars.seed << '\n';
            ui->output->setPlainText(QString::fromStdString(s.str()));
        }

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
            QApplication::processEvents();
            std::stringstream s;
            s << "t: " << t;
            ui->output->appendPlainText(QString::fromStdString(s.str()));

            if (t == 0) metapop.exitburnin();

            // Life cycle of the metapopulation
            metapop.cycle(pars, arch);

            // Is the population still there?
            if (metapop.isextinct()) {
                std::cout << "The population went extinct at t = " << t << '\n';
                break;
            }

            // Analyze the metapopulation if needed
            //if (timetosave(t, pars)) {
                collector.analyze(metapop, pars, arch);

                std::vector<double> fst_vals = collector.get_Fst();
                plot_fst(fst_vals);


                std::vector<double> gst_vals = collector.get_Gst();
                plot_gst(gst_vals);

                std::vector<double> eco_trait_vals = collector.get_eco_trait(metapop);
                std::vector<double> sex_trait_vals = collector.get_sex_trait(metapop);
                std::vector<double> neu_trait_vals = collector.get_neu_trait(metapop);

                plot_hist(ui->plot_eco_trait, ecoBars, eco_trait_vals);
                plot_hist(ui->plot_sex_trait, sexBars, sex_trait_vals);
                plot_hist(ui->plot_neu_trait, neuBars, neu_trait_vals);
        }
        // Show output
        {
            std::stringstream s;
            s
              << "EI: " << collector.getEI() << '\n'
              << "RI: " << collector.getRI() << '\n'
              << "SI: " << collector.getSI() << '\n'
            ;
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
    is_running = false;
}
