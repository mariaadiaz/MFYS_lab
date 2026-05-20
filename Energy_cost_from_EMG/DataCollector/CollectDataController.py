"""
Controller class for the Data Collector GUI
This is the controller for the GUI that lets you connect to a base, scan via rf for sensors, and stream data from them in real time.
"""
# Maybe useful link: https://stackoverflow.com/questions/13302393/qobjectconnect-cannot-queue-arguments-of-type-object-in-pyside

import pickle
import time  # For timestamping
from collections import deque
import threading
import pandas as pd
from itertools import zip_longest
from Functions.synergies import *
import matplotlib.pyplot as plt
from datetime import datetime
from PyQt5.QtWidgets import QMessageBox, QDialog, QTableWidgetItem

from Plotter.GenericPlot import *
from Plotter.PlotOutput import *
from AeroPy.TrignoBase import *
from AeroPy.DataManager import *

clr.AddReference("System.Collections")

base = TrignoBase()
TrigBase = base.BaseInstance

app.use_app('PyQt5')


class PlottingManagement:
    def __init__(self, metrics, emgplot=None, process=False, processingpanel = None):
        self.process = process
        self.EMGplot = emgplot
        self.metrics = metrics
        self.processpanel  = processingpanel # Panel to update the variables from the new interface
        self.ProcessingHandler = Processingcollector() # Functions used to calculate the synergies and power outputs
        self.packetCount = 0  # Number of packets received from base
        self.pauseFlag = True  # Flag to start/stop collection and plotting
        self.numSamples = 10000  # Default number of samples to visualize at a time
        self.DataHandler = DataKernel(TrigBase)  # Data handler for receiving data from base
        self.outData = [[0]]
        self.Index = None
        self.newTransform = None

        # Timestamp path and filename to save
        timestr = time.strftime("%d%m%Y_%H%M")
        # self.path_to_save = timestr + "_participantX.pkl"
        # self.file_name = timestr + "_participantX_times.txt"

        # Define the path to save the collected data and the file name to save the starting and stopping times. For testing purposes, I am using fixed names but in the future I will use the timestr variable to have unique files for each session.
        self.path_to_save = 'Testing_recent.pkl'
        self.file_name = "Testing.txt"
        self.min2 = True
        self.min3 = True
        self.iterations = 1
        self.counter = 1
        self.sample_x = np.array([75, 95, 126])
        self.sample_y = np.array([])  # Initialize for both live and uploaded data modes

        # Initialize the max values of the muscles to zero. Afterwards,the peaks of the muscles will be compared to these max values to update them and to normalize the muscle activations for the synergy extraction and power calculations
        self.max_rta, self.max_lta, self.max_rrf, self.max_lrf = 0,0,0,0
        self.max_rgm, self.max_lgm, self.max_rgl, self.max_lgl = 0,0,0,0

        # ############ Only for testing with previously collected data ############
        # with open('Data\\Participant 22\\08122023_P22_02_EMG.pkl', 'rb') as f:
        #     self.data = pickle.load(f)
        # self.data.remove(self.data[32])
        # self.values = 3 * 60 * np.round(1925.9259259259259)
        # self.start = 0
        # #########################################################################

    def streaming(self):
        """This is the data processing thread"""
        self.emg_queue = deque()
        while self.pauseFlag is True:
            continue
        while self.pauseFlag is False:
            self.DataHandler.processData(self.emg_plot)
            self.updatemetrics()

    def analysis_data(self):
        """
        Processing thread — triggers at 80s (min2) and 140s (min3) each cycle.

        Two data modes (switch by commenting/uncommenting):
          1) Previously collected data  — uses self.data loaded from pickle
          2) Real-time data             — uses self.DataHandler.getStoredData()
        """        
        while self.process:

            if self.DataHandler.restart_var == True:
                self.min2 = True
                self.min3 = True
                self.DataHandler.restart_var = False
                
            if self.DataHandler.timer > 80 and self.min2 == True:
                self.min2 = False
                ################################
                # 1) Previously collected data #
                ################################

                # # t_emg = np.linspace(0,len(self.data)/np.round(self.sampleRates[0]),len(self.data))
                # self.fin = int(self.values * (self.iterations - 1) + self.values)
                # partial_emg = [vals[self.start:self.fin] for vals in self.data] # 3 minutes every time
                # self.init_emg = [vals[0:int(self.values * 3)] for vals in self.data]
                # x = [vals[0:int((1.334) * 60 * np.round(self.sampleRates[0]))] for vals in partial_emg] # 20 seconds every time
                # y = pd.DataFrame(list(zip_longest(*x)))

                ############################
                # 2) Real-time data
                # (comment out block above and uncomment below)
                ############################
                [x,n_a] = self.DataHandler.getStoredData()
                flat_data = self.make_data_flat(x)
                y = pd.DataFrame(list(zip_longest(*flat_data)))
                ############################
                
                y.columns = self.channelnames
                [power1, W1, H1, VAF, self.max_rta, self.max_lta, self.max_rrf, self.max_lrf, self.max_rgm, self.max_lgm, self.max_rgl, self.max_lgl] = self.ProcessingHandler.process_emg(y, np.round(self.sampleRates[0]), self.max_rta, self.max_lta, self.max_rrf, self.max_lrf, self.max_rgm, self.max_lgm, self.max_rgl, self.max_lgl)
                print("min 2 results (power & VAF):", np.round(np.sum(power1),3), np.round(VAF,3))
                

            if self.DataHandler.timer > 140 and self.min3 == True:            
                self.min3 = False

                ################################
                # 1) Previously collected data #
                ################################
                # x = [vals[int((1.334) * 60 * np.round(self.sampleRates[0])):int((2.334) * 60 * np.round(self.sampleRates[0]))] for vals in partial_emg]
                # # x = [vals[0:int((self.counter + 1.334) * 60 * np.round(self.sampleRates[0]))] for vals in partial_emg]
                # y = pd.DataFrame(list(zip_longest(*x)))
                # self.start = self.fin
                ############################
                # 2) Real-time data
                # (comment out block above and uncomment below)
                ############################
                [n_a, x] = self.DataHandler.getStoredData()
                flat_data = self.make_data_flat(x)
                y = pd.DataFrame(list(zip_longest(*flat_data)))
                ############################

                y.columns = self.channelnames
                [power2, W2, H2, VAF2, self.max_rta, self.max_lta, self.max_rrf, self.max_lrf, self.max_rgm, self.max_lgm, self.max_rgl, self.max_lgl] = self.ProcessingHandler.process_emg(y, np.round(self.sampleRates[0]), self.max_rta, self.max_lta, self.max_rrf, self.max_lrf, self.max_rgm, self.max_lgm, self.max_rgl, self.max_lgl)
                [mean_power, mean_ssv, mean_dist, EC_estimation, n_H2] = self.ProcessingHandler.compare_muscle_outputs(power1, power2, W1, W2, H1, H2, n=4)

                print("min 3 results (power & VAF):", np.round(np.sum(power2),3), np.round(VAF2,3))
                
                self.sample_y = np.append(self.sample_y,EC_estimation)
                print("this is self_sample_y:", self.sample_y)
                print("EC estimation:", np.round(EC_estimation,3))
                print("Iteration {}".format(self.iterations),": Mean power, Mean ssv, Mean dist:", np.round(mean_power,2), np.round(mean_ssv,2), np.round(mean_dist,2))

                """Update the variables and the canvas"""
                self.updatevariables(np.round(mean_power,2),np.round(mean_ssv,2),np.round(mean_dist,2),np.round(EC_estimation,3))           
                self.updatecanvas(H1, n_H2)


                # Optimization section
                if self.iterations >= 3:
                    if self.DataHandler.initialization == True:
                        self.DataHandler.initialization = False

                    [y_pred, y_std, ei, new_x, x_range] = self.ProcessingHandler.bayesian_optimization(self.sample_x, self.sample_y)
                    if self.iterations == 3:
                        first_pred = y_pred
                    
                    val = np.where(x_range == self.closest_value(x_range, new_x))
                    new_y = first_pred[val]

                    self.updateoptimization(self.iterations, x_range, y_pred, y_std, ei, self.sample_x, self.sample_y, new_x, new_y[0])

                    # This condition is to check if the optimization has converged
                    if 0.95 * self.sample_x[-1]  <= np.round(new_x) <= 1.05 * self.sample_x[-1] :
                        # Calculate the mean between new_x and the previous sf value
                        mean_sf = np.mean([self.sample_x[-1], new_x])
                        print("The optimization has converged you can stop")
                        # print("The final EC value is:", np.round(self.sample_y[-1]))
                        # last element from sample_y is the EC value
                        print("The optimal SF is:", np.round(mean_sf))
                    if any(0.99 * value <= np.round(new_x) <= 1.01 * value for value in self.sample_x):
                        print("Stop: new SF is in the 1% range of some values in sample_x")
                        print("The optimal SF is:", np.round(new_x))
                    self.sample_x = np.append(self.sample_x, np.round(new_x))
                    print("SF values:", self.sample_x)

                self.iterations += 1
                # self.counter += 3
    
    def closest_value(self, vector, x):
        closest_index = np.argmin(np.abs(vector - x))
        return vector[closest_index]                  

    def vispyPlot(self):
        """Real-time EMG plot thread."""
        skipCount = 0
        while self.pauseFlag is False:
            if len(self.emg_plot) >= 2:
                incData = self.emg_plot.popleft()  # Data at time T-1
                self.outData = list(np.asarray(incData, dtype='object')[tuple([self.dataStreamIdx])])
                if self.dataStreamIdx and self.outData[0].size > 0:
                    try:
                        self.EMGplot.plot_new_data(self.outData, [self.emg_plot[0][i][0] for i in self.dataStreamIdx])
                    except IndexError:
                        print("Index Error Occurred: CollectDataController.py - line 102")

    def make_data_flat(self, x):
        """Flatten allcollectiondata list for real-time processing path."""
        process_data = list()
        for obj in x:
            process_data.append(obj)

        partial_data = list(tuple(zip(*process_data)))

        flat_data = list()
        for ii in range(len(partial_data)):
            flat_data.append([item for sublist in partial_data[ii] for item in sublist])
        return flat_data

    def updatemetrics(self):
        self.metrics.framescollected.setText(str(self.DataHandler.packetCount))
    
    def updatevariables(self,var1,var2,var3,var4):
        # Results panel
        self.processpanel.power_output.setText(str(var1))
        self.processpanel.similarity_output.setText(str(var2))
        self.processpanel.distance_output.setText(str(var3))
        # EC estimation
        self.processpanel.ec_output.setText(str(var4))
    
    def updatecanvas(self, H1, H2):
        # print("Update canvas was reached")
        self.processpanel.update_plot(H1, H2)

    def updateoptimization(self, iter, x_range, y_pred, y_std, ei, sample_x, sample_y, new_x, new_y):
        # Next SF value
        self.processpanel.nextvalue_output.setText(str(np.round(new_x)))
        # Update optimization panel
        self.processpanel.update_optimization_plot(iter, x_range, y_pred, y_std, ei, sample_x, sample_y, new_x, new_y)

    def threadManager(self, start_trigger, stop_trigger):
        """Handles the threads for the DataCollector gui"""
        self.emg_plot = deque()

        self.t1 = threading.Thread(target=self.streaming)
        self.t1.setDaemon(True)
        self.t1.start()

        if self.EMGplot:
            self.t2 = threading.Thread(target=self.vispyPlot)
            self.t2.setDaemon(True)
            if not start_trigger:
                self.t2.start()

        if start_trigger:
            self.t3 = threading.Thread(target=self.WaitingForStartTrigger)
            self.t3.start()

        if stop_trigger:
            self.t4 = threading.Thread(target=self.WaitingForStopTrigger)
            self.t4.start()
        
        if self.process:
            self.t5 = threading.Thread(target=self.analysis_data)
            self.t5.start()


    def WaitingForStartTrigger(self):
        while TrigBase.IsWaitingForStartTrigger():
            continue
        self.pauseFlag = False
        if self.EMGplot:
            self.t2.start()
        print("Trigger Start - Collection Started")

    def WaitingForStopTrigger(self):
        while TrigBase.IsWaitingForStartTrigger():
            continue
        while TrigBase.IsWaitingForStopTrigger():
            continue
        self.pauseFlag = True
        self.metrics.pipelinestatelabel.setText(self.PipelineState_Callback())
        print("Trigger Stop - Data Collection Complete")
        self.DataHandler.processData(self.emg_plot)

    # ---------------------------------------------------------------------------------
    # ---- Callback Functions
    def PipelineState_Callback(self):
        return TrigBase.GetPipelineState()

    def Connect_Callback(self):
        """Callback to connect to the base"""
        TrigBase.ValidateBase(key, license)

    def FrameCount_Callback(self):
        return self.DataHandler.packetCount

    def Pair_Callback(self):
        """Callback to tell the base to enter pair mode for new sensors"""
        TrigBase.PairSensor()
        dialogbox = QDialog()
        dialogbox.setWindowTitle("Awaiting sensor pair request. . .")
        dialogbox.setFixedWidth(300)
        dialogbox.setFixedHeight(50)
        dialogbox.show()
        while TrigBase.CheckPairStatus():
            continue
        dialogbox.close()
        messagebox = QMessageBox()
        messagebox.setText("Pair another sensor?")
        messagebox.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        messagebox.setIcon(QMessageBox.Question)
        button = messagebox.exec_()

        if button == QMessageBox.Yes:
            self.Pair_Callback()
        else:
            self.Scan_Callback()

    def Scan_Callback(self):
        """Callback to tell the base to scan for any available sensors"""
        f = TrigBase.ScanSensors().Result
        self.nameList = TrigBase.GetSensorNames()
        self.SensorsFound = len(self.nameList)

        # TrigBase.SelectSensor(0)
        TrigBase.SelectAllSensors()
        return self.nameList

    def Start_Callback(self, start_trigger, stop_trigger):
        """Callback to start the data stream from Sensors"""
        column = 0
        if not start_trigger:
            self.pauseFlag = False
        self.metrics.framescollected.setText("0")
        self.DataHandler.packetCount = 0
        self.DataHandler.allcollectiondata = [[]]

        if TrigBase.GetPipelineState() == 'Armed':
            for i in range(len(self.channelnames)):
                self.DataHandler.allcollectiondata.append([])

        if TrigBase.GetPipelineState() == 'Connected':
            self.channelcount = 0
            TrigBase.Configure(start_trigger, stop_trigger)

            self.channelnames = []
            self.sampleRates = []
            self.samplesPerFrame = [[] for i in range(self.SensorsFound)]
            self.variables_names = [[] for i in range(self.SensorsFound)]
            self.plotCount = 0
            # ---- Discover sensor channels
            self.dataStreamIdx = []  # This list indexes into the sensor data array, selecting relevant data to visualize
            idxVal = 0

            for i in range(self.SensorsFound):
                selectedSensor = TrigBase.GetSensorObject(i)
                if len(selectedSensor.TrignoChannels) > 0:
                    for channel in range(len(selectedSensor.TrignoChannels)):
                        self.channelcount += 1
                        self.DataHandler.allcollectiondata.append([])
                        self.channelnames.append(selectedSensor.TrignoChannels[channel].Name)
                        self.sampleRates.append(selectedSensor.TrignoChannels[channel].SampleRate)
                        self.samplesPerFrame.append(selectedSensor.TrignoChannels[channel].SamplesPerFrame)
                        # ---- Collect the EMG channels for visualization, excluding skin check channels
                        if "EMG" in selectedSensor.TrignoChannels[channel].Name:
                            if "TrignoAvanti" in str(selectedSensor) and "2" in selectedSensor.TrignoChannels[
                                channel].Name:  # Avanti skin check
                                pass
                            elif "AvantiDoubleMini" in str(selectedSensor) and "2" in selectedSensor.TrignoChannels[
                                channel].Name:  # Duo Mini skin check 1
                                pass
                            elif "AvantiDoubleMini" in str(selectedSensor) and "4" in selectedSensor.TrignoChannels[
                                channel].Name:  # Duo Mini skin check 2
                                pass
                            else:
                                self.dataStreamIdx.append(idxVal)
                                self.plotCount += 1

                        if self.channelnames[idxVal] == 'EMG 1':
                            if selectedSensor.PairNumber == 1 or selectedSensor.PairNumber == 9:
                                self.channelnames[idxVal] = 'RightRectusFemoris EMG'

                            elif selectedSensor.PairNumber == 2 or selectedSensor.PairNumber == 10:
                                self.channelnames[idxVal] = 'RightGastroLateralis EMG'

                            elif selectedSensor.PairNumber == 3 or selectedSensor.PairNumber == 11:
                                self.channelnames[idxVal] = 'RightGastroMedialis EMG'

                            elif selectedSensor.PairNumber == 4 or selectedSensor.PairNumber == 12:
                                self.channelnames[idxVal] = 'RightTibialisAnterior EMG'

                            elif selectedSensor.PairNumber == 5 or selectedSensor.PairNumber == 13:
                                self.channelnames[idxVal] = 'LeftRectusFemoris EMG'

                            elif selectedSensor.PairNumber == 6 or selectedSensor.PairNumber == 14:
                                self.channelnames[idxVal] = 'LeftGastroLateralis EMG'

                            elif selectedSensor.PairNumber == 7 or selectedSensor.PairNumber == 15:
                                self.channelnames[idxVal] = 'LeftGastroMedialis EMG'

                            elif selectedSensor.PairNumber == 8 or selectedSensor.PairNumber == 16:
                                self.channelnames[idxVal] = 'LeftTibialisAnterior EMG'
                            else:
                                pass
                        idxVal += 1
        self.metrics.totalchannels.setText(str(self.channelcount))
        if self.EMGplot:
            self.EMGplot.initiateCanvas(None, None, self.plotCount, 1, self.numSamples)
        TrigBase.Start()
        self.threadManager(start_trigger, stop_trigger)

        ########## Starting time ##########
        # Print and save the starting time
        now = datetime.now()
        print("Starting time =", now)
        with open(self.file_name, "w") as file:
            file.write("Starting time:" + str(now))

    def Stop_Callback(self):
        """Callback to stop the data stream"""
        self.pauseFlag = True
        TrigBase.Stop()
        print("Data Collection Complete")
        datos = self.DataHandler.allcollectiondata
        with open(self.path_to_save, 'wb') as f:
             pickle.dump(datos, f)
        
        """" To load the file in another ocasion I have to use the lines below"""
        # with open('all_data.pkl', 'rb') as f:
        #     mynewlist = pickle.load(f)

        ########## Stopping time ##########
        # Print and save the starting time
        now = datetime.now()
        print("Stopping time =", now)
        with open(self.file_name, "a") as file:
            file.write("\nStopping time:" + str(now))

    # ---------------------------------------------------------------------------------
    # ---- Helper Functions
    def getSampleModes(self, sensorIdx):
        """Gets the list of sample modes available for selected sensor"""
        sampleModes = TrigBase.AvailibleSensorModes(sensorIdx)
        return sampleModes

    def getCurMode(self, sensorIdx):
        """Gets the current mode of the sensors"""
        curModes = TrigBase.GetCurrentSensorMode(sensorIdx)
        return curModes

    def setSampleMode(self, curSensor, setMode):
        """Sets the sample mode for the selected sensor"""
        TrigBase.SetSampleMode(curSensor, setMode)