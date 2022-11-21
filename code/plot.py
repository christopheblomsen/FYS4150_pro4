import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import os

# from matplotlib.ticker import ScalarFormatter
plt.rc('text', usetex=True)
# plt.rc('font', family='Helvetica')

SMALL_SIZE = 15
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# plt.style.use('seaborn')
# plt.rcParams.update({'font.size': 30})

# params = {
#     'legend.fontsize': 9,
#     'font.size': 9,
#     'figure.figsize': (8.647/2, 8.0/2),
#     'axes.labelsize': 9,
#     'axes.titlesize': 9,
#     'xtick.labelsize': 9,
#     'ytick.labelsize': 9
#     }

# plt.rcParams.update(params)

def compute_heat_capacity(energy, T, L, kb=1):
    return 1 / (L*L) * 1 / (kb * T*T) * (np.var(energy))

def compute_m_abs(magnetisation, L):
    return np.mean(np.abs(magnetisation)) / (L*L)

def compute_epsilon(energy, L):
    return np.mean(energy) / (L*L)

def compute_susceptibility(magnetisation, T, L, kb=1):
    return 1 / (L*L) * 1 / (kb * T) * np.var(np.abs(magnetisation))


def Z(beta):
    return 12 + 4*np.cosh(8*beta)

def avg_E(beta):
    return -32*np.sinh(8*beta) / Z(beta)

def avg_eps(beta):
    return avg_E(beta)/4 

def avg_E2(beta):
    return 256/Z(beta) * np.cosh(8*beta) 

def avg_M(beta):
    return 8/Z(beta) * (2 + np.exp(8*beta))

def avg_M2(beta):
    return 32/Z(beta) * (1 + np.exp(8*beta))

def CV(beta):
    return 1/4 * (avg_E2(beta) - avg_E(beta)**2)

def chi(beta):
    return 1/4 * (avg_M2(beta) - avg_M(beta)**2)

beta = 1 



os.chdir(r"/home/rhuvy/Documents/FYS4150/Projects/FYS4150_pro4/code")

# Create the arrays for problem 4
epsilon = np.zeros(20)
m = np.zeros(20)
cv = np.zeros(20)
Chi = np.zeros(20)

# Create the arrays for problem 5
epsilon1_ord = np.zeros(20)
m1_ord = np.zeros(20)

epsilon1_rdm = np.zeros(20)
m1_rdm = np.zeros(20)

epsilon2_ord = np.zeros(20)
m2_ord = np.zeros(20)

epsilon2_rdm = np.zeros(20)
m2_rdm = np.zeros(20)

# list for every number of cycle computed
list = ["500", "1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000"
, "10000", "20000", "30000", "40000", "50000", "60000", "70000", "80000", "90000", "100000"]

# to access the correct location in the array
index = 0

for i in list:
    # Problem 4
    # data = pa.mat()    
    # data.load(f"L2_T1_{i}_cycles.bin")
    
    # epsilon[index] = compute_epsilon(data[0,:], 2)
    # m[index] = compute_m_abs(data[1,:], 2)
    # cv[index] = compute_heat_capacity(data[0,:], 1., 2)
    # Chi[index] = compute_susceptibility(data[1,:], 1., 2)

    # Problem 5
    data1_ord = pa.mat()
    data1_rdm = pa.mat()
    data2_rdm = pa.mat()
    data2_ord = pa.mat()

    data1_ord.load(f"L20_ord_T1_{i}_cycles.bin")
    data1_rdm.load(f"L20_random_T1_{i}_cycles.bin")
    data2_ord.load(f"L20_ord_T24_{i}_cycles.bin")
    data2_rdm.load(f"L20_random_T24_{i}_cycles.bin")

    epsilon1_ord[index] = compute_epsilon(data1_ord[0,:], 20)
    epsilon1_rdm[index] = compute_epsilon(data1_rdm[0,:], 20)
    epsilon2_ord[index] = compute_epsilon(data2_ord[0,:], 20)
    epsilon2_rdm[index] = compute_epsilon(data2_rdm[0,:], 20)

    m1_ord[index] = compute_m_abs(data1_ord[1,:], 20)
    m1_rdm[index] = compute_m_abs(data1_rdm[1,:], 20)
    m2_ord[index] = compute_m_abs(data2_ord[1,:], 20)
    m2_rdm[index] = compute_m_abs(data2_rdm[1,:], 20)
    

    index += 1


# Plot Problem 4
test = list[::3]

# fig, ax = plt.subplots()
# ax.plot(list, Chi, 'ro-', label='cv')
# ax.hlines(chi(1)/4, list[0], list[-1], 'k', linestyles='dashed', label='cv exact')
# ax.set_xticks(test)
# ax.set_xticklabels(test)
# ax.legend(loc='best')
# ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
# ax.set_ylabel(r'$\mathrm{Susceptibility}$')
# ax.set_title(r'$\mathrm{Susceptibility\ vs\ number\ of\ cycles}$')
# fig.savefig("chi_problem4.pdf")

# fig, ax = plt.subplots()
# ax.plot(list, cv, 'ro-', label='cv')
# ax.hlines(CV(1), list[0], list[-1], 'k', linestyles='dashed', label='cv exact')
# ax.set_xticks(test)
# ax.set_xticklabels(test)
# ax.legend(loc='best')
# ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
# ax.set_ylabel(r'$\mathrm{C_v}$')
# ax.set_title(r'$\mathrm{C_v\ vs\ number\ of\ cycles}$')
# fig.savefig("CV_problem4.pdf")
# plt.show()

# Plot Problem 5

fig, ax = plt.subplots()
ax.plot(list, epsilon1_ord, 'bo-', label=r'$\mathrm{Ordered, T=1}$')
ax.plot(list, epsilon1_rdm, 'go-', label=r'$\mathrm{Random, T=1}$')
ax.plot(list, epsilon2_ord, 'b--', label=r'$\mathrm{Ordered, T=2.4}$')
ax.plot(list, epsilon2_rdm, 'g--', label=r'$\mathrm{Random, T=2.4}$')
ax.set_xticks(test)
ax.set_xticklabels(test)
ax.legend(loc='best')
ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
ax.set_ylabel(r'$<\epsilon>$')
ax.set_title(r'$<\epsilon> \mathrm{\ vs\ number\ of\ cycles}$')
fig.savefig("epsilon_problem5.pdf")

fig, ax = plt.subplots()
ax.plot(list, epsilon1_ord, 'bo-', label=r'$\mathrm{Ordered, T=1}$')
ax.plot(list, epsilon1_rdm, 'go-', label=r'$\mathrm{Random, T=1}$')
ax.plot(list, epsilon2_ord, 'b--', label=r'$\mathrm{Ordered, T=2.4}$')
ax.plot(list, epsilon2_rdm, 'g--', label=r'$\mathrm{Random, T=2.4}$')
ax.set_xticks(test)
ax.set_xticklabels(test)
ax.legend(loc='best')
ax.set_xlabel(r'$<|m|>\mathrm{Number\ of\ cycles}$')
ax.set_ylabel(r'$<|m|>$')
ax.set_title(r'$<|m|> \mathrm{\ vs\ number\ of\ cycles}$')
fig.savefig("m_problem5.pdf")


plt.show()

# data1 = pa.mat()
# data1.load("test_serial_0.bin")
# data1 = np.array(data1)

# data2 = pa.mat()
# data2.load("test_serial_1.bin")
# data2 = np.array(data2)

# data3 = pa.mat()
# data3.load("test_serial_2.bin")
# data3 = np.array(data3)

# data4 = pa.mat()
# data4.load("test_serial_3.bin")
# data4 = np.array(data4)

# data5 = pa.mat()
# data5.load("test_serial_4.bin")
# data5 = np.array(data5)

# data6 = pa.mat()
# data6.load("test_serial_5.bin")
# data6 = np.array(data6)

# data7 = pa.mat()
# data7.load("test_serial_6.bin")
# data7 = np.array(data7)

# data8 = pa.mat()
# data8.load("test_serial_7.bin")
# data8 = np.array(data8)

# data9 = pa.mat()
# data9.load("test_serial_8.bin")
# data9 = np.array(data9)

# data10 = pa.mat()
# data10.load("test_serial_9.bin")
# data10 = np.array(data10)

# data11 = pa.mat()
# data11.load("test_serial_10.bin")
# data11 = np.array(data11)

# data12 = pa.mat()
# data12.load("test_serial_11.bin")
# data12 = np.array(data12)

# data13 = pa.mat()
# data13.load("test_serial_12.bin")
# data13 = np.array(data13)

# data14 = pa.mat()
# data14.load("test_serial_13.bin")
# data14 = np.array(data14)

# data15 = pa.mat()
# data15.load("test_serial_14.bin")
# data15 = np.array(data15)

# data16 = pa.mat()
# data16.load("test_serial_15.bin")
# data16 = np.array(data16)

# data17 = pa.mat()
# data17.load("test_serial_16.bin")
# data17 = np.array(data17)

# data18 = pa.mat()
# data18.load("test_serial_17.bin")
# data18 = np.array(data18)

# data19 = pa.mat()
# data19.load("test_serial_18.bin")
# data19 = np.array(data19)

# data20 = pa.mat()
# data20.load("test_serial_19.bin")
# data20 = np.array(data20)