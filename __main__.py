from Model.data import Data
from Model.analysis import Process

print("initializing analysis class")
data = Data()
proc = Process()

data.read_data()

proc.powerspectrum(data)
