"""
Data processing GUI
This is the GUI that displays the synergies and power outputs from the EMG signals.
"""
import sys
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
# from Functions.synergies import *
import tkinter as tk
from Plotter import PlotOutput as gp

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class ProcessingWindow(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.dataprocessingpanel = self.processingpanel()


    def processingpanel(self):
        # QWidget.__init__(self)
        # panel = QWidget()
        layout = QHBoxLayout()
        self.setStyleSheet("background-color:lavender;")

        self.power_val = ''
        self.similarity_val = ''
        self.distance_val =''
        # self.fullresultspanel = self.completeresults()
        self.nextvalue = self.relevant_outcomes()
        self.plotPanel = self.Plotter()
        self.resultspanel = self.Lastresultspanel()
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.plotPanel)
        self.splitter.addWidget(self.resultspanel)

        self.splitter2 = QSplitter(Qt.Horizontal)
        # self.splitter2.addWidget(self.fullresultspanel)
        self.splitter2.addWidget(self.nextvalue)
        self.splitter2.addWidget(self.splitter)


        layout.addWidget(self.splitter2)
        self.setLayout(layout)
        self.setWindowTitle("Process Data GUI")
        # processingdatapanel.setFixedWidth(150)

        return self
    
    #-----------------------------------------------------------------------
    #---- GUI Components
    def relevant_outcomes(self):

        valuepanel = QWidget()
        valuepanel_layout = QVBoxLayout()

        self.next_value = QLabel('Next SF Value')
        self.next_value.setAlignment(Qt.AlignCenter)
        self.next_value.setStyleSheet("border :2px solid black;background : cornflowerblue")
        valuepanel_layout.addWidget(self.next_value)

        # Increase font size and make it bold
        font = self.next_value.font()
        font.setPointSize(20)  # Set your desired font size
        font.setBold(True)
        self.next_value.setFont(font)


        self.nextvalue_output = QLabel('-')
        self.nextvalue_output.setAlignment(Qt.AlignCenter)
        self.nextvalue_output.setStyleSheet("color:black")
        valuepanel_layout.addWidget(self.nextvalue_output)

        # Increase font size and make it bold
        font_output = self.nextvalue_output.font()
        font_output.setPointSize(20)  # Set your desired font size
        font_output.setBold(True)
        self.nextvalue_output.setFont(font_output)

        self.ec_estimation = QLabel('EC Estimation')
        self.ec_estimation.setAlignment(Qt.AlignCenter)
        self.ec_estimation.setStyleSheet("border :2px solid black;background : mediumturquoise")
        valuepanel_layout.addWidget(self.ec_estimation)
        
        # Increase font size and make it bold
        font = self.ec_estimation.font()
        font.setPointSize(20)  # Set your desired font size
        font.setBold(True)
        self.ec_estimation.setFont(font)

        self.ec_output = QLabel('-')
        self.ec_output.setAlignment(Qt.AlignCenter)
        self.ec_output.setStyleSheet("color:black")
        valuepanel_layout.addWidget(self.ec_output)

        # Increase font size and make it bold
        font_output = self.ec_output.font()
        font_output.setPointSize(20)  # Set your desired font size
        font_output.setBold(True)
        self.ec_output.setFont(font_output)

        valuepanel.setLayout(valuepanel_layout)
        return valuepanel
    
    def completeresults(self):
        self.tableWidget = QTableWidget()

        #Row count
        rows = 10
        self.tableWidget.setRowCount(rows)

        #Column count
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setHorizontalHeaderItem(0, QTableWidgetItem("Power"))
        self.tableWidget.setHorizontalHeaderItem(1, QTableWidgetItem("Similarity"))
        self.tableWidget.setHorizontalHeaderItem(2, QTableWidgetItem("Distance"))
        
        for i in range(rows):
            self.tableWidget.setItem(i,0,QTableWidgetItem('-'))
            self.tableWidget.setItem(i,1,QTableWidgetItem('-'))
            self.tableWidget.setItem(i,2,QTableWidgetItem('-'))
        
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
        self.power.setStyleSheet("border :2px solid black;background : lightyellow")
        power_layout.addWidget(self.power)

        self.power_output = QLabel('-')
        self.power_output.setAlignment(Qt.AlignCenter)
        self.power_output.setStyleSheet("color:black")
        power_layout.addWidget(self.power_output)
        resultspanel_layout.addLayout(power_layout)

       # Increase font size for power output
        font_power_label = self.power.font()
        font_power_label.setPointSize(12)  # Set your desired font size
        font_power_label.setBold(True)
        self.power.setFont(font_power_label)
        font_power = self.power_output.font()
        font_power.setPointSize(12)  # Set your desired font size
        font_power.setBold(True)
        self.power_output.setFont(font_power)

        # similarity
        self.similarity = QLabel('Similarity')
        self.similarity.setAlignment(Qt.AlignCenter)
        self.similarity.setStyleSheet("border :2px solid black;background:lightblue")
        similarity_layout.addWidget(self.similarity)

        self.similarity_output = QLabel('-')
        self.similarity_output.setAlignment(Qt.AlignCenter)
        self.similarity_output.setStyleSheet("color:black")
        similarity_layout.addWidget(self.similarity_output)
        resultspanel_layout.addLayout(similarity_layout)

       # Increase font size for similarity output and label
        font_similarity_label = self.similarity.font()
        font_similarity_label.setPointSize(12)  # Set your desired font size
        font_similarity_label.setBold(True)
        self.similarity.setFont(font_similarity_label)
        font_similarity = self.similarity_output.font()
        font_similarity.setPointSize(12)  # Set your desired font size
        font_similarity.setBold(True)
        self.similarity_output.setFont(font_similarity)

        # Distance
        self.distance = QLabel('Distance')
        self.distance.setAlignment(Qt.AlignCenter)
        self.distance.setStyleSheet("border :2px solid black;background:lightcoral")
        distance_layout.addWidget(self.distance)

        self.distance_output = QLabel('-')
        self.distance_output.setAlignment(Qt.AlignCenter)
        self.distance_output.setStyleSheet("color:black")
        distance_layout.addWidget(self.distance_output)
        resultspanel_layout.addLayout(distance_layout)

       # Increase font size for distance output and label
        font_distance_label = self.distance.font()
        font_distance_label.setPointSize(12)  # Set your desired font size
        font_distance_label.setBold(True)
        self.distance.setFont(font_distance_label)
        font_distance = self.distance_output.font()
        font_distance.setPointSize(12)  # Set your desired font size
        font_distance.setBold(True)
        self.distance_output.setFont(font_distance)

        resultspanel.setLayout(resultspanel_layout)

        return resultspanel

    def Plotter(self):
        widget = QWidget()
        widget.setLayout(QVBoxLayout())
        self.canvas_wrapper = gp.CanvasWrapper()
        widget.layout().addWidget(self.canvas_wrapper.canvas.native)
        return widget
    
    def update_plot(self,x,x2):
        self.canvas_wrapper.updating_canvas(x,x2)

    def update_optimization_plot(self, iter, x_range, y_pred, y_std, ei, sample_x, sample_y, new_x, new_y):
        self.canvas_wrapper.optimization_canvas(iter, x_range, y_pred, y_std, ei, sample_x, sample_y, new_x, new_y)


