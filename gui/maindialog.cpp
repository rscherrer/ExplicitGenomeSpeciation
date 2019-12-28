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
            if (timetosave(t, pars)) collector.analyze(metapop, pars, arch);

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
