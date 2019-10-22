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
  ui->plot->xAxis->setRange(0,90);
  ui->plot->yAxis->setRange(0,1);
  ui->plot->xAxis->grid()->setVisible(false);
  ui->plot->yAxis->grid()->setVisible(false);
 // plot->xAxis->grid().setVisible(false);
}

MainDialog::~MainDialog()
{
    delete ui;
}


void MainDialog::plot()
{
    ui->plot->graph(0)->clearData();
    ui->plot->clearItems();
    ui->plot->graph(0)->setData(qv_x, qv_y);
    ui->plot->replot();
    ui->plot->update();
//    QThread::msleep(100);
   // ui->plot->
}

Param MainDialog::createPars()
{
    Param pars;
    pars.seed = static_cast<size_t>(ui->rng_seed->value());
    return pars;
}

void MainDialog::on_run_button_clicked()
{
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
        //for (int t = -pars.tburnin; t < pars.tend; ++t) {

        for(int t = 0; t < 10; ++t) {
            std::stringstream s;
            s << "t: " << t;
            ui->output->appendPlainText(QString::fromStdString(s.str()));
            //if (t == 0) metapop.exitburnin();

            // Life cycle of the metapopulation
            //metapop.cycle(pars, arch);

            // Is the population still there?
            if (metapop.isextinct()) {
                std::cout << "The population went extinct at t = " << t << '\n';
                break;
            }

            // Analyze the metapopulation if needed
           // if (timetosave(t, pars)) {
               // collector.analyze(metapop, pars, arch);
                //std::vector<double> fst_vals = collector.get_Fst();
                qv_x.clear();
                qv_y.clear();
                for(size_t i = 0; i < 90; ++i) {
                    qv_x.append(i);
                    qv_y.append(t);
                }
               // plot();

              //  ui->plot->clearGraphs();
              //  ui->plot->addGraph();
               // ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
              //  ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone);
             //   ui->plot->graph(0)->clearData();

                ui->plot->graph(0)->setData(qv_x, qv_y);
              //  ui->plot->graph(0)->rescaleAxes();
                ui->plot->replot(QCustomPlot::rpQueued);
                ui->plot->update();
                QThread::msleep(100);
            // }
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

void MainDialog::on_btn_rand_points_clicked()
{
    for(int t = 0; t < 1000; ++t) {
        qv_x.clear();
        qv_y.clear();
        for(size_t i = 0; i < 90; ++i) {
            qv_x.append(i);
            qv_y.append(1.0 * std::rand() / RAND_MAX);
        }

      //  ui->plot->graph(0)->clearData();
      //  ui->plot->clearItems();
      //  ui->plot->graph(0)->setData(qv_x, qv_y);

       // ui->plot->replot();
     //   ui->plot->update();




        std::stringstream s;
        s << t << "\n";
        ui->output->appendPlainText(QString::fromStdString(s.str()));
        ui->output->repaint();
        ui->output->update();
        QThread::msleep(10);
    }
}



void MainDialog::on_btn_one_point_clicked()
{
    double x = 90.0 * std::rand() / RAND_MAX;
    double y = 1.0 * std::rand() / RAND_MAX;
    //qv_x.clear();
    qv_x.append(x);
    //qv_y.clear();
    qv_y.append(y);

    //ui->plot->graph(0)->clearData();
   // ui->plot->clearItems();
    ui->plot->graph(0)->setData(qv_x, qv_y);
    ui->plot->replot();
    ui->plot->update();
}
