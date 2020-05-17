#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include <QDialog>
#include "gui/qcustomplot.h"
#include "library/Param.h"

namespace Ui {
  class MainDialog;
}

class MainDialog : public QDialog
{
  Q_OBJECT

public:
  explicit MainDialog(QWidget *parent = nullptr);
  ~MainDialog();

  void plot_fst(const std::vector<double>& v);
  void plot_gst(const std::vector<double>& v);
  void plot_eco_trait(const std::vector<double>& v);

  void update_plot_popsize(int t, size_t n_0, size_t n_1);



private slots:
  void on_run_button_clicked();

  void on_pushButton_clicked();

private:
  Ui::MainDialog *ui;

  ///Read the parameters from the GUI
  Param createPars();
  void setup_spinboxes();

  QVector<double> pop_x;
  QVector<double> pop_y_0;
  QVector<double> pop_y_1;

  bool is_running;

  QCPBars *ecoBars_0;
  QCPBars *ecoBars_1;

  QCPBars *sexBars_0;
  QCPBars *sexBars_1;

  QCPBars *neuBars_0;
  QCPBars *neuBars_1;
};

#endif // MAINDIALOG_H
