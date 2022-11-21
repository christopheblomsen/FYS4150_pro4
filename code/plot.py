import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

data = pa.mat()

data.load("test_serial1.000000.bin")
data = np.array(data)
test = data[1,:]
print(np.mean(test))
