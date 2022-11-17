import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

data = pa.mat()

N = np.linspace(0, 2, 2)
data.load('output.bin')
print(data[0, :])

eps = data[0, :]
mag = data[1, :]
N = data[2, :]

plt.plot(N, eps)
plt.xlabel('N')
plt.ylabel(r'$\epsilon$')
plt.show()

plt.plot(N, mag)
plt.xlabel('N')
plt.ylabel(r'$\epsilon$')
plt.show()
