import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

data = pa.mat()

N = np.linspace(0, 2, 2)
data.load('output.bin')
print(data[0, :])

plt.plot(data[2, :], data[0, :])
plt.show()
