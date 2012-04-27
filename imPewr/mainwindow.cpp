//#include "H5Cpp.h"
#include "mainwindow.h"
#include "imagepane.h"
#include "Array.h"
#include <vector>


using namespace Array;

// Specify what 
// typedef array3<Real>    StackReal;


MainWindow::MainWindow()
	: QMainWindow()
{
	// Set up window appearance
	setMinimumSize(800,600);

	// Set up upper toolbar
	//addToolBar(Qt::ToolBarArea Qt::TopToolBarArea,topMenuBar);
	fileMenu = menuBar()->addMenu("File");
	// Open dataset
	openDatasetAction = new QAction("Open Dataset",fileMenu);
	fileMenu->addAction(openDatasetAction);
	//connect(openDatasetAction,SIGNAL(triggered()),this,SLOT
	// Open image series
	openImageSeriesAction = new QAction("Open Series of Images",fileMenu);
	fileMenu->addAction(openImageSeriesAction);
	//connect
	// Save dataset
	saveDatasetAction = new QAction("Save Dataset",fileMenu);
	fileMenu->addAction(saveDatasetAction);
	// Close whole program
	closeAction = new QAction("Exit Program",fileMenu);
	fileMenu->addAction(closeAction);
	connect(closeAction,SIGNAL(triggered()),this,SLOT(close()));


	// Set up image controls pane
	QDockWidget *imageControlsWidget = new QDockWidget(this);
	QLabel *temp1 = new QLabel("Wheee!",imageControlsWidget);
	addDockWidget(Qt::LeftDockWidgetArea,imageControlsWidget);


	// Set up image pane
	ImagePane *imagePane = new ImagePane(this);
	QLabel *temp2 = new QLabel("Wheee2!",imagePane);
	QLabel *imageSliderBar = new QLabel("SliderBar");
	QVBoxLayout *imagePaneLayout = new QVBoxLayout();
	QWidget *imageWidget = new QWidget();
	imagePaneLayout->addWidget(imagePane);
	imagePaneLayout->addWidget(imageSliderBar);
	imageWidget->setLayout(imagePaneLayout);
	setCentralWidget(imageWidget);


	// Set up properties pane
	QDockWidget *propertiesWidget = new QDockWidget(this);
	QLabel *temp3 = new QLabel("Wheee3!",propertiesWidget);
	addDockWidget(Qt::RightDockWidgetArea,propertiesWidget);



	// Set up all actions that can be performed on images or datasets
	createActions();



	// Set up main window layout
	//QHBoxLayout *mainLayout = new QHBoxLayout();
	//mainLayout->addWidget(imageControlsWidget);
	//mainLayout->addWidget(imagePane);
	//mainLayout->addWidget(propertiesWidget);
	//setLayout(mainLayout);


};


void MainWindow::createActions()
{
	openImageSeriesAction = new QAction(tr("&Open..."), this);
	openImageSeriesAction->setShortcut(tr("Ctrl+O"));
	connect(openImageSeriesAction, SIGNAL(triggered()), this, SLOT(openImageSeries()));
};


void MainWindow::openImageSeries()
{
};


struct MicrographStack {
	int				xlength, ylength, zlength;
	array3<float>	data;
};



//struct micrographStack {
//	
//}
	//	array3<float>	data;
//
//}


