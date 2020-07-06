import os
import numpy as np
from nptdms import TdmsFile as tdms

try:
    from tkinter import Tk
    from tkFileDialog import askopenfilenames
except:
    from tkinter import Tk
    from tkinter import filedialog

#####################################################

class Data:

    def __init__(self):
        self.prop = {}
        self.data = {}

    def read_data(self):
        Tk().withdraw()

        self.prop['fns'] = filedialog.askopenfilenames(filetypes=[("TDMS","*.tdms")])
        self.prop['Nf'] = len(self.prop['fns'])

        mdata = self.load_mdata()
        self.check_files(mdata)
        print("file check ok")

        self.prop['samples'] = list(mdata[0].items())[0][1]
        self.prop['rate'] = list(mdata[0].items())[1][1]
        self.prop['Tmsr'] = self.prop['samples']/self.prop['rate']
        self.prop['Nchannels'] = len(list(mdata[0].items())[3][1])


        for i in range(self.prop['Nf']):
            data = tdms.read(self.prop['fns'][i])
            for j in range(self.prop['Nchannels']):
                if i == 0:
                    self.data[j] = np.array(data.groups()[0].channels()[j].data)
                else:
                    self.data[j] = np.c_[self.data[j], np.array(data.groups()[0].channels()[j].data)]


    def load_mdata(self):
        mdata = {}
        for i in range(self.prop['Nf']):
            temp = tdms.read_metadata(self.prop['fns'][i])
            mdata[i] = temp.groups()[0].properties
        return mdata

    def load_data(self):
        for i in range(self.prop['Nf']):
            if i == 0:
                data = tdms.read(self.prop['fns'][i])
            else:
                data = np.c_[data, tdms.read(self.prop['fns'][i])]
        return data

    def check_files(self, mdata):
        for i in range(self.prop['Nf']):
            if i == self.prop['Nf']-1:
                break
            else:
                if (mdata[i] == mdata[i+1]):
                    continue
                else:
                    print('mdata from multiple selected files do not correspond')
                exit()

    def split_traces(self, Nsplit):
        if (self.prop['Nf'] > 1):
            print('Error: multiple time traces selected, or time trace already split')
            return

        if ((self.prop['Tmsr']/Nsplit).is_integer() == False):
            print('Error: choos Nsplit for which your splitted Tmsr is an integer number')
        else:
            Chunks = int(self.prop['samples']/Nsplit)
            for j in range(self.prop['Nchannels']):
                temp = np.empty([0, Chunks])
                for i in range(0, self.prop['samples'], Chunks):
                    if i == 0:
                        temp = self.data[j][i:i+Chunks]
                    else:
                        temp = np.c_[temp, self.data[j][i:i+Chunks]]
                self.data[j] = temp

        self.prop['Nf'] = Nsplit
        self.prop['Tmsr'] = self.prop['Tmsr']/Nsplit
