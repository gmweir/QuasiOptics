import os, sys, time
import numpy as np
import datetime

"""
Description:
class progressbar provides a customizable progress indicator

Author:
Fabian Wilde, IPP HGW

Date:
September 2016
"""

class progressBar(object):
    def __init__(self, **kwargs):
        if 'char' in kwargs.keys():
            self.char = kwargs['char']
        else:
            self.char = '#'
        if 'barlen' in kwargs.keys():
            if not type(kwargs['len']) is int:
                raise ValueError(" expecting integer for input argument 'len'.")
            self.barlen = kwargs['len']
        else:
            self.barlen = 50
        if 'bar' in kwargs.keys():
            if not kwargs['bar'] in [True, False]:
                raise ValueError(" expecting boolean or string representing boolean keyword for input argument 'bar'.")
            self.bar = kwargs['bar']
        else:
            self.bar = True
        if 'eta' in kwargs.keys():
            if not kwargs['eta'] in [True, False]:
                raise ValueError(" expecting boolean or string representing boolean keyword for input argument 'eta'.")
            self.eta = kwargs['eta']
        else:
            self.eta = True
        if 'spinner' in kwargs.keys():
            if not kwargs['spinner'] in [True, False]:
                raise ValueError(" expecting boolean or string representing boolean keyword for input argument 'spinner'.")
            self.spinner = kwargs['spinner']
        else:
            self.spinner = True
        if 'spinnerdat' in kwargs.keys():
            self.spinnerdat = kwargs['spinnerdat']
        else:
            self.spinnerdat = ['|', '/', '-', '\\']
        if 'value' in kwargs.keys():
            if not type(kwargs['value']) is int:
                raise ValueError(" expecting integer for input argument 'value'.")
            self.value = kwargs['value']
        else:
            self.value = 0
        if 'total' in kwargs.keys():
            if not type(kwargs['total']) is int:
                raise ValueError(" expecting integer for input argument 'total'.")
            self.total = kwargs['total']
        else:
            self.total = 100
        if 'text' in kwargs.keys():
            self.customText = kwargs['text']
        else:
            self.customText = 'Progress:'

        self.progress = round((float(self.value) / float(self.total)) * 100, 2)
        self.timeLeft = -1
        self.times = []
        self.values = []
        self.spinnerCount = 0

    def show(self):
        str2show = self.customText + ' ' + str(self.progress) + " % "
        if self.bar:
            # build progressbar
            numchars = int((self.progress / 100) * self.barlen)
            str2show += "["
            # add bar characters according to progress value
            for i in range(0, numchars):
                str2show = str2show + self.char
            # fill up with white space
            if (self.barlen - numchars) > 0:
                for i in range(0, self.barlen - numchars):
                    str2show += " "
            str2show += "]"

        if self.eta:
            if self.timeLeft == -1:
                str2show += " - ETA: TBD"
            elif np.isnan(self.timeLeft):
                str2show += " - ETA: TBD"
            else:
                hh = int(self.timeLeft / 3600)
                mm = int((self.timeLeft - hh * 3600) / 60)
                ss = int((self.timeLeft - hh * 3600 - mm * 60))
                str2show += " - ETA: " + str(hh) + ":" + str(mm) + ":" + str(ss)

        if self.spinner:
            self.spinnerCount += 1
            if self.spinnerCount > (len(self.spinnerdat) - 1):
                self.spinnerCount = 0
            str2show += " " + self.spinnerdat[self.spinnerCount] + " "

        print(str2show, end='\r\n')

    def update(self, *kwargs):
        value = kwargs[0]
        if len(kwargs) > 1:
            self.ramUsage = kwargs[1]
        self.value = value
        self.values.append(value)
        self.times.append(time.time())
        if len(self.values) > 1:
            uniqueValues=np.unique(self.values)
            uniqueTimes = []
            for index in range(0,len(uniqueValues)):
                uniqueTimes.append(self.times[np.where(np.array(self.values) == uniqueValues[index])[0][0]])
            self.values = list(uniqueValues)
            self.times = list(uniqueTimes)
        if self.value > self.total:
            self.value = self.total
        if len(self.times) >= 2:
            self.timeLeft = (self.total - self.value)/(np.mean(np.diff(self.values))/np.mean(np.diff(self.times)))
        else:
            self.timeLeft = -1
        tmp = round((float(self.value) / float(self.total)) * 100, 2)
        if tmp >= 0:
            self.progress = tmp
        self.show()

    def finish(self,timeElapsed=False):
        self.update(self.total)
        print("")
        if timeElapsed:
            print(str(self.times[-1]-self.times[0])+" seconds elapsed.")

    def clear(self):
        self.times = []
        self.timeLeft = -1
        self.progress = 0
