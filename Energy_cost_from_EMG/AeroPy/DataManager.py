"""
This is the class that handles the data that is output from the Delsys Trigno Base.
Create an instance of this and pass it a reference to the Trigno base for initialization.
See CollectDataController.py for a usage example.
"""

from collections import deque
import numpy as np
import time

class DataKernel():
    def __init__(self, trigno_base):
        self.TrigBase = trigno_base
        self.packetCount = 0
        self.sampleCount = 0
        self.allcollectiondata = [[]]
        self.count_exception = 0
        self.init = True
        self.thedata = deque()
        self.data_min2 = deque()
        self.data_min3 = deque()
        self.start_time = time.perf_counter()
        self.restart_var = False
        self.initialization = True
        self.data_initialization = deque()

    def processData(self, data_queue):
        if self.init:
            self.init = False
            self.start_time = time.perf_counter()
            
        self.timer = time.perf_counter() - self.start_time
        """Processes the data from the Trigno Base and places it in the data_queue argument"""
        outArr = self.GetData()
        if outArr is not None:
            for i in range(len(outArr)):
                self.allcollectiondata[i].extend(outArr[i][0].tolist())
            try:
                for i in range(len(outArr[0])):
                    if np.asarray(outArr[0]).ndim == 1:
                        data_queue.append(list(np.asarray(outArr, dtype='object')[0]))
                        # self.thedata.append(list(np.asarray(outArr, dtype='object')[0]))
                        if (self.timer <= 80): # 20 seconds after min 1
                            self.data_min2.append(list(np.asarray(outArr, dtype='object')[0]))
                        if (self.timer >= 80 and self.timer <= 140): # From min 2 to min 3
                            self.data_min3.append(list(np.asarray(outArr, dtype='object')[0]))
                        if self.initialization == True:
                            self.data_initialization.append(list(np.asarray(outArr, dtype='object')[0]))

                    else:
                        data_queue.append(list(np.asarray(outArr, dtype='object')[:, i]))
                        # self.thedata.append(list(np.asarray(outArr, dtype='object')[:, i]))
                        if (self.timer <= 80): # 20 seconds after min 1
                            self.data_min2.append(list(np.asarray(outArr, dtype='object')[:,i]))
                        if (self.timer > 80 and self.timer <= 140): # From min 2 to min 3
                            self.data_min3.append(list(np.asarray(outArr, dtype='object')[:,i]))
                        if self.initialization == True:
                            self.data_initialization.append(list(np.asarray(outArr, dtype='object')[:,i]))

                
                '''Re-start the timer'''
                if self.timer >= 180:
                    print("Timer was restarted")
                    self.start_time = time.perf_counter()
                    self.restart_var = True
                    # Re-initialize
                    self.data_min2 = deque()
                    self.data_min3 = deque()

                try:
                    self.packetCount += len(outArr[0])
                    self.sampleCount += len(outArr[0][0])
                except:
                    pass
            except IndexError:
                pass

    def GetData(self):
        """Dictionary: Callback to get the data from the streaming sensors"""
        dataReady = self.TrigBase.CheckDataQueue()
        if dataReady:
            try:
                DataOut = self.TrigBase.PollData()
                if len(DataOut) > 0:
                    outArr = [[] for i in range(len(DataOut.Keys))]
                    keys = []
                    for i in DataOut.Keys:
                        keys.append(i)
                    for j in range(len(DataOut.Keys)):
                        outBuf = DataOut[keys[j]]
                        outArr[j].append(np.asarray(outBuf, dtype='object'))
                    return outArr
                else:
                    return None
            except:
                print("Exception in DataOut")
                self.count_exception =+ 1
                pass
        else:
            return None


    # -----------------------------------------------------------
    # ---- Helper Functions
    def getPacketCount(self):
        return self.packetCount

    def resetPacketCount(self):
        self.packetCount = 0
        self.sampleCount = 0

    def getSampleCount(self):
        return self.sampleCount

    # ------------------------------------------------------------
    # ---- Functions to store the data
    def getStoredData(self):
        # return self.thedata
        return self.data_min2, self.data_min3
