from DataCollector.CollectDataWindow import CollectDataWindow
from DataCollector.ProcessingWindow_v2 import ProcessingWindow
from StartMenu.StartWindow import StartWindow
import sys
from PyQt5.QtWidgets import *

class FrameController():
    def __init__(self):
        self.startWindow = StartWindow(self)
        self.collectWindow = CollectDataWindow(self)
        self.processWindow = ProcessingWindow()

        self.startWindow.show()

        self.curHeight = 650
        self.curWidth = 1115

    def showStartMenu(self):
        self.collectWindow.close()
        self.startWindow.show()

    def showCollectData(self):
        self.startWindow.close()
        if self.startWindow.plot_enabled.isChecked():
            self.collectWindow.plot_enabled = True
            self.collectWindow.AddPlotPanel()
        if self.startWindow.process_enabled.isChecked():
            self.collectWindow.process_enabled = True
            # self.processWindow.show()
        self.collectWindow.SetCallbackConnector()
        self.collectWindow.show()


def main():
    app = QApplication(sys.argv)
    controller = FrameController()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()




