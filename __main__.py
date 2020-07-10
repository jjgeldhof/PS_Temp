from Model.data import Data
from Model.analysis import Process
import numpy as np

print("initializing analysis class")
dat = Data()
self = Process()

dat.read_data()

self.powerspectrum(dat)

self.bin_spectrum(900)
self.fit()
