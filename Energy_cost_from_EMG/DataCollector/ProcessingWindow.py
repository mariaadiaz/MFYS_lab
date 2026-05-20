"""
Data processing GUI
This is the GUI that displays the synergies and power outputs from the EMG signals.
"""
import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
# from Functions.synergies import *
from DataCollector.CollectDataController import *
import tkinter as tk
from Plotter import GenericPlot as gp

class ProcessingWindow(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        # self.processingpanel = self.dataprocessingpanel()

    # def dataprocessingpanel(self):
        # QWidget.__init__(self)
        self.power_val = ''
        self.similarity_val = ''
        self.distance_val =''
        self.fullresultspanel = self.completeresults()
        self.plotPanel = self.Plotter()
        self.resultspanel = self.Lastresultspanel()
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.plotPanel)
        self.splitter.addWidget(self.resultspanel)

        self.splitter2 = QSplitter(Qt.Horizontal)
        self.splitter2.addWidget(self.fullresultspanel)
        self.splitter2.addWidget(self.splitter)

        layout = QHBoxLayout()
        self.setStyleSheet("background-color:#3d4c51;")
        layout.addWidget(self.splitter2)
        self.setLayout(layout)
        self.setWindowTitle("Process Data GUI")

        # self.valx1 = 'No value yet'
        # self.valx2 = 'No value yet'
        # self.valx3 = 'No value yet'        
        # self.timer = QTimer(self, interval = 30 * 1000) # every 30*1000 ms
        # self.timer.timeout.connect(self.panelupdatevalues)
        # self.timer.start()
    
        # return self
        #---- Connect the controller to the GUI
        # self.CallbackConnector = PlottingManagement(None, None, None)

    # def panelupdatevalues(self):
    #     print('updating values')
    #     self.power_output.setText(str(self.valx1))
    #     self.similarity_output.setText(str(self.valx2))
    #     self.distance_output.setText(str(self.valx3))
        
        # self.getpowervalue(self.valx1)
        # self.getsimilarityvalue(self.valx2)
        # self.getdistancevalue(self.valx3)        
    #-----------------------------------------------------------------------
    #---- GUI Components
    def completeresults(self):
        self.tableWidget = QTableWidget()

        #Row count
        self.tableWidget.setRowCount(10)

        #Column count
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setHorizontalHeaderItem(0, QTableWidgetItem("Power"))
        self.tableWidget.setHorizontalHeaderItem(1, QTableWidgetItem("Similarity"))
        self.tableWidget.setHorizontalHeaderItem(2, QTableWidgetItem("Distance"))
        # self.tableWidget.setItem(1,1,10)
        
        return self.tableWidget

    def Lastresultspanel(self):

        resultspanel = QWidget()
        resultspanel_layout = QHBoxLayout()

        power_layout = QVBoxLayout()
        similarity_layout = QVBoxLayout()
        distance_layout = QVBoxLayout()
        # Power
        self.power = QLabel('Power')
        self.power.setAlignment(Qt.AlignCenter)
        self.power.setStyleSheet("border :2px solid black;background : yellow")
 
        power_layout.addWidget(self.power)

        self.power_output = QLabel('-')
        self.power_output.setAlignment(Qt.AlignCenter)
        self.power_output.setStyleSheet("color:white")
        power_layout.addWidget(self.power_output)
        resultspanel_layout.addLayout(power_layout)

        # similarity
        self.similarity = QLabel('Similarity')
        self.similarity.setAlignment(Qt.AlignCenter)
        self.similarity.setStyleSheet("border :2px solid black;background : blue")
        similarity_layout.addWidget(self.similarity)

        self.similarity_output = QLabel('-')
        self.similarity_output.setAlignment(Qt.AlignCenter)
        self.similarity_output.setStyleSheet("color:white")
        similarity_layout.addWidget(self.similarity_output)
        resultspanel_layout.addLayout(similarity_layout)

        # Distance
        self.distance = QLabel('Distance')
        self.distance.setAlignment(Qt.AlignCenter)
        self.distance.setStyleSheet("border :2px solid black;background : red")
        distance_layout.addWidget(self.distance)

        self.distance_output = QLabel('-')
        self.distance_output.setAlignment(Qt.AlignCenter)
        self.distance_output.setStyleSheet("color:white")
        distance_layout.addWidget(self.distance_output)
        resultspanel_layout.addLayout(distance_layout)

        resultspanel.setLayout(resultspanel_layout)

        return resultspanel

    def Plotter(self):
        widget = QWidget()
        widget.setLayout(QVBoxLayout())
        plot_mode = 'windowed'                 # Select between 'scrolling' and 'windowed'
        pc = gp.GenericPlot(plot_mode)
        pc.native.objectName = 'vispyCanvas'
        pc.native.parent = self
        widget.layout().addWidget(pc.native)
        self.plotCanvas = pc
        return widget

    #-----------------------------------------------------------------------
    #---- Callback Functions
    def getpowervalue(self,x):
        # self.power_val = self.CallbackConnector.Partial_results_callback()
        self.power_output.setText(str(x))

    def getsimilarityvalue(self,x):
        # self.similarity_val = self.CallbackConnector.Partial_results_callback()
        self.similarity_output.setText(str(x))

    def getdistancevalue(self,x):
        # self.distance_val = self.CallbackConnector.Partial_results_callback()
        self.distance_output.setText(str(x))

 

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ProcessingWindow = ProcessingWindow()
    ProcessingWindow.show()
    sys.exit(app.exec_())