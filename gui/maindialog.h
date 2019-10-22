#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include <QDialog>
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

  void plot();

private slots:
  void on_run_button_clicked();

  void on_btn_rand_points_clicked();

  void on_btn_one_point_clicked();

private:
  Ui::MainDialog *ui;

  ///Read the parameters from the GUI
  Param createPars();

  QVector<double> qv_x;
  QVector<double> qv_y;
};

#endif // MAINDIALOG_H
