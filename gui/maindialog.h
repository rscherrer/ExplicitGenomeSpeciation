#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include <QDialog>
#include "library/Param.h"
#include "qcustomplot.h"

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

  void update_plot_popsize(int t, size_t n, size_t n_0, size_t n_1);



private slots:
  void on_run_button_clicked();

  void on_pushButton_clicked();

private:
  Ui::MainDialog *ui;

  ///Read the parameters from the GUI
  Param createPars();
  void setup_spinboxes();

  QVector<double> pop_x;
  QVector<double> pop_y;
  QVector<double> pop_y_0;
  QVector<double> pop_y_1;

  bool is_running;

  QCPBars *ecoBars;
  QCPBars *sexBars;
  QCPBars *neuBars;
};

#endif // MAINDIALOG_H
