#include "maindialog.h"
#include "ui_maindialog.h"
#include "library/Simul.h"
#include <sstream>

MainDialog::MainDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::MainDialog)
{
  ui->setupUi(this);
}

MainDialog::~MainDialog()
{
  delete ui;
}

Param MainDialog::createPars()
{
    Param pars;
    pars.seed = static_cast<size_t>(ui->rng_seed->value());
    return pars;
}

void set_pars_thijs(Param& pars) {
    pars.rdynamics = 1;
    pars.inflow = 1;
    pars.outflow = 1000;
    pars.capacity = 1000;
    pars.replenish = 1.0;
    pars.hsymmetry = 0;
    pars.ecosel = 1;
    pars.dispersal = 0.001;
    pars.birth = 4;
    pars.survival = 0.6;
    pars.sexsel = 10;
    pars.matingcost = 0.01;
    pars.maxfeed = 0.0004;
    return;
}



void MainDialog::on_run_button_clicked()
{
    try
    {
        // Create the parameters from the GUI
        Param pars = createPars();
        set_pars_thijs(pars);

        QScatterSeries *series0 = new QScatterSeries();
        series0->setName("FST");
        series0->setMarkerSize(5.0);

        QChart *chart = new QChart();
        chart->addSeries(series0);

        QValueAxis *axisX = new QValueAxis();
        axisX->setRange(0, 300);

        QValueAxis *axisY = new QValueAxis();
        axisY->setRange(-0.01, 0.01);


        chart->addAxis(axisX, Qt::AlignBottom);
        chart->addAxis(axisY, Qt::AlignLeft);
        series0->attachAxis(axisY);
        series0->attachAxis(axisX);
        chart->setDropShadowEnabled(false);

        ui->graphicsView->setChart(chart);


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
                collector.analyze(metapop, pars);

                std::vector<float> to_plot = collector.get_Fst();
                // fun update
        //        for(size_t i = 0; i < to_plot.size(); ++i) {
          //          to_plot[i] = 1.0 * std::rand() / RAND_MAX;
           //     }


                // visualize output
                series0->clear();
                for(size_t i = 0; i < to_plot.size(); ++i) {
                    *series0 << QPointF( i, to_plot[i]);
                }



                ui->graphicsView->repaint();
            //  }
                std::stringstream s;
                s << "t: " << t;
                ui->output->appendPlainText(QString::fromStdString(s.str()));


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

          //  ui->graphicsView->chart()
        }
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
    }
}
