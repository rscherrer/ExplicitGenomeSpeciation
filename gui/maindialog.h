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

private slots:
  void on_run_button_clicked();

private:
  Ui::MainDialog *ui;

  ///Read the parameters from the GUI
  Param createPars();
};

#endif // MAINDIALOG_H
