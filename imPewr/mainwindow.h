#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>
//#include <QMainWindow>


class QAction;


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

	

private slots:
	openImageSeries();
	//QFileDialog dialog(this);
	//dialog.setFileMode(QFileDialog::AnyFile);


private:
	// Top toolbar
	// File menu
	QMenu *fileMenu;
	QAction *openDatasetAction;
	QAction *openImageSeriesAction;
	QAction *saveDatasetAction;
	QAction *closeAction;
	//QLabel *testing;


	// Image data
	// vector<ArrayReal *> imageStack;
	//QFileDialog *filename;
	void createActions();
	QAction *openImageSeriesAction;
	


signals:


//private slots:
	

};



#endif