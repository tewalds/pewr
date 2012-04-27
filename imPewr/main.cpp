#include <QApplication>
#include "mainwindow.h"


 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
	 
     MainWindow MainWindow;
	 MainWindow.show();


//	 QPushButton *button = new QPushButton(
//         QApplication::translate("childwidget", "Press me"), &MainWindow);
//     button->move(100, 100);
//     button->show();

     return app.exec();
 }