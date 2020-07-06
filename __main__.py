from Model.data import Data
from Model.analysis import Process

print("initializing analysis class")
dat = Data()
proc = Process()

dat.read_data()

proc.powerspectrum(dat)
