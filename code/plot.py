import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import os
from scipy.stats import linregress

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
    return np.var(energy) / (L*L * kb * T*T) 

def compute_m_abs(magnetisation, L):
    return np.mean(np.abs(magnetisation)) / (L*L) 

def compute_epsilon(energy, L):
    return np.mean(energy) / (L*L)

def compute_susceptibility(magnetisation, T, L, kb=1):
    # m2 = np.mean(magnetisation**2)
    # m = np.mean(np.abs(magnetisation))
    return np.var(np.abs(magnetisation)) / (L*L * kb * T)





def Z(beta):
    return 12 + 4*np.cosh(8*beta)

def avg_E(beta):
    return -32*np.sinh(8*beta) / Z(beta)

def avg_eps(beta):
    return avg_E(beta)/4 

def avg_E2(beta):
    return 256/Z(beta) * np.cosh(8*beta) 

def avg_M(beta):
    return 32/Z(beta) * (2 + np.exp(8*beta))

def avg_M2(beta):
    return 32/Z(beta) * (1 + np.exp(8*beta))

def CV(beta):
    return 1/4 * (avg_E2(beta) - avg_E(beta)**2)

def chi(beta):
    return 1/4 * (avg_M2(beta) - avg_M(beta)**2)

beta = 1 



#os.chdir(r"/home/rhuvy/Documents/FYS4150/Projects/FYS4150_pro4/code/data")
os.chdir(r"data")

# # Create the arrays for problem 4
# epsilon = np.zeros(20)
# m = np.zeros(20)
# cv = np.zeros(20) 
# Chi = np.zeros(20)

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
list_ = ["500", "1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000"
, "10000", "20000", "30000", "40000", "50000", "60000", "70000", "80000", "90000", "100000"]

# # to access the correct location in the array
index = 0

for i in list_:
#     # Problem 4
#     data = pa.mat()    
#     data.load(f"L2_T1_{i}_cycles.bin")
    
#     epsilon[index] = compute_epsilon(data[0,:], 2)
#     m[index] = compute_m_abs(data[1,:], 2)
#     cv[index] = compute_heat_capacity(data[0,:], 1., 2)
#     Chi[index] = compute_susceptibility(data[1,:], 1., 2)

#     # Problem 5
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


# # Plot Problem 4
test = list_[::3]

# fig, ax = plt.subplots()
# ax.plot(list_, Chi, 'ro-', label='$\chi$')
# ax.hlines(chi(1)/4, list_[0], list_[-1], 'k', linestyles='dashed', label=r'$\chi$ exact')
# ax.set_xticks(test)
# ax.set_xticklabels(test)
# ax.legend(loc='best')
# ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
# ax.set_ylabel(r'$\mathrm{Susceptibility}$')
# ax.set_title(r'$\mathrm{Susceptibility\ vs\ number\ of\ cycles}$')
# # fig.savefig("chi_problem4.pdf")

# fig, ax = plt.subplots()
# ax.plot(list_, cv, 'ro-', label='cv')
# ax.hlines(CV(1), list_[0], list_[-1], 'k', linestyles='dashed', label='cv exact')
# ax.set_xticks(test)
# ax.set_xticklabels(test)
# ax.legend(loc='best')
# ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
# ax.set_ylabel(r'$\mathrm{C_v}$')
# ax.set_title(r'$\mathrm{C_v\ vs\ number\ of\ cycles}$')
# # fig.savefig("CV_problem4.pdf")
# plt.show()

# Plot Problem 5

# fig, ax = plt.subplots()
# ax.plot(list_, epsilon1_ord, 'bo-', label=r'$\mathrm{Ordered, T=1}$')
# ax.plot(list_, epsilon1_rdm, 'go-', label=r'$\mathrm{Random, T=1}$')
# ax.plot(list_, epsilon2_ord, 'b--', label=r'$\mathrm{Ordered, T=2.4}$')
# ax.plot(list_, epsilon2_rdm, 'g--', label=r'$\mathrm{Random, T=2.4}$')
# ax.set_xticks(test)
# ax.set_xticklabels(test)
# ax.legend(loc='best')
# ax.set_xlabel(r'$\mathrm{Number\ of\ cycles}$')
# ax.set_ylabel(r'$<\epsilon>$')
# ax.set_title(r'$<\epsilon> \mathrm{\ vs\ number\ of\ cycles}$')
# # fig.savefig("epsilon_problem5.pdf")

fig, ax = plt.subplots()
ax.plot(list_, m1_ord, 'bo-', label=r'$\mathrm{Ordered, T=1}$')
ax.plot(list_, m1_rdm, 'go-', label=r'$\mathrm{Random, T=1}$')
ax.plot(list_, m2_ord, 'b--', label=r'$\mathrm{Ordered, T=2.4}$')
ax.plot(list_, m2_rdm, 'g--', label=r'$\mathrm{Random, T=2.4}$')
ax.set_xticks(test)
ax.set_xticklabels(test)
ax.legend(loc='best')
ax.set_xlabel(r'$<|m|>\mathrm{Number\ of\ cycles}$')
ax.set_ylabel(r'$<|m|>$')
ax.set_title(r'$<|m|> \mathrm{\ vs\ number\ of\ cycles}$')
# # fig.savefig("m_problem5.pdf")


# Problem 6
# Normally after the iteration the latest saved data1_rdm and data2_rdm 
# will be L=20, N=100000

# eps_prob6 = data1_rdm[0,:] /400 
# eps2_prob6 = data2_rdm[0,:] / 400

# var_prob5 = np.round(np.var(eps_prob6), 5)
# var2_prob5 = np.round(np.var(eps2_prob6), 5)

# plt.figure()
# plt.hist(eps_prob6, density=True, bins=[-2,-1.99,-1.98,-1.96, -1.94], label=f"Var: {var_prob5}")
# plt.xlabel(r"$\epsilon$")
# plt.ylabel("Density")
# plt.legend(loc='best')
# plt.title(r"Normalized histogram of $\epsilon$ for T=1")
# # plt.savefig("Histo_T1.pdf")

# plt.figure()
# plt.hist(eps2_prob6, density=True, bins=20, label=f"Var: {var2_prob5}")
# plt.xlabel(r"$\epsilon$")
# plt.ylabel("Density")
# plt.legend(loc='best')
# plt.title(r"Normalized histogram of $\epsilon$ for T=2.4")
# # plt.savefig("Histo_T2.pdf")


# # Problem 8 and 9 
# temperature = pa.mat()
# temperature.load("temp_prob8.bin")


# taille = 10

# temperature = np.array(temperature).reshape(taille)

# eps40 = np.zeros(taille)
# m40 = np.zeros(taille)
# cv40 = np.zeros(taille)
# chi40 = np.zeros(taille)

# eps60 = np.zeros(taille)
# m60 = np.zeros(taille)
# cv60 = np.zeros(taille)
# chi60 = np.zeros(taille)

# eps80 = np.zeros(taille)
# m80 = np.zeros(taille)
# cv80 = np.zeros(taille)
# chi80 = np.zeros(taille)

# eps100 = np.zeros(taille)
# m100 = np.zeros(taille)
# cv100 = np.zeros(taille)
# chi100 = np.zeros(taille)

# for k in range(len(temperature)):
#     T = temperature[k]
#     # Load the data for one temperature
#     data8_40 = pa.mat()
#     data8_40.load(f"L2_prob8_{k}.bin")
#     data8_40 = np.array(data8_40)

#     # data8_60 = pa.mat()
#     # data8_60.load(f"L60_prob8_{k}.bin")
#     # data8_60 = np.array(data8_60)

#     # data8_80 = pa.mat()
#     # data8_80.load(f"L80_prob8_{k}.bin")
#     # data8_80 = np.array(data8_80)
    
#     # data8_100 = pa.mat()
#     # data8_100.load(f"L100_prob8_{k}.bin")
#     # data8_100 = np.array(data8_100)

#     # Calculate the four quantities
#     eps40[k] = compute_epsilon(data8_40[0,:], 2.)
#     # m40[k] = compute_m_abs(data8_40[1,:], 40.)
#     # cv40[k] = compute_heat_capacity(data8_40[0,:],  T, 40.)
#     # chi40[k] = compute_susceptibility(data8_40[1,:],  T, 40.)

#     # eps60[k] = compute_epsilon(data8_60[0,:], 60.)
#     # m60[k] = compute_m_abs(data8_60[1,:], 60.)
#     # cv60[k] = compute_heat_capacity(data8_60[0,:],  T, 60.)
#     # chi60[k] = compute_susceptibility(data8_60[1,:],  T, 60.)

#     # eps80[k] = compute_epsilon(data8_80[0,:], 80.)
#     # m80[k] = compute_m_abs(data8_80[1,:], 80.)
#     # cv80[k] = compute_heat_capacity(data8_80[0,:],  T, 80.)
#     # chi80[k] = compute_susceptibility(data8_80[1,:],  T, 80.)

#     # eps100[k] = compute_epsilon(data8_100[0,:], 100.)
#     # m100[k] = compute_m_abs(data8_100[1,:], 100.)
#     # cv100[k] = compute_heat_capacity(data8_100[0,:],  T, 100.)
#     # chi100[k] = compute_susceptibility(data8_100[1,:],  T, 100.)

# # import sys
# # np.set_printoptions(threshold=sys.maxsize)

# # print(data8_40[0,:])
# # plt.figure()
# # plt.plot(temperature, cv100, 'k', label="100")
# # plt.plot(temperature, cv80, 'o', label="80")
# # plt.plot(temperature, cv60, 'b', label="60")
# # plt.plot(temperature, cv40, 'g', label="40")
# # plt.legend(loc='best')
# # plt.xlabel("Temperature [J/k]")
# # plt.ylabel(r"$C_v$")

# # plt.figure()
# # plt.plot(temperature, chi100, 'k', label="100")
# # plt.plot(temperature, chi80, 'o', label="80")
# # plt.plot(temperature, chi60, 'b', label="60")
# # plt.plot(temperature, chi40, 'g', label="40")
# # plt.legend(loc='best')
# # plt.xlabel("Temperature [J/k]")
# # plt.ylabel(r"$\chi$")

# plt.figure()
# # plt.plot(temperature, m100, 'k', label="100")
# # plt.plot(temperature, m80, 'o', label="80")
# # plt.plot(temperature, m60, 'b', label="60")
# plt.plot(temperature, m40, 'g', label="40")
# plt.legend(loc='best')
# plt.xlabel("Temperature [J/k]")
# plt.ylabel(r"$|m|$")

# # plt.figure()
# # plt.plot(temperature, eps100, 'k', label="100")
# # plt.plot(temperature, eps80, 'o', label="80")
# # plt.plot(temperature, eps60, 'b', label="60")
# # plt.plot(temperature, eps40, 'g', label="40")
# # plt.legend(loc='best')
# # plt.xlabel("Temperature [J/k]")
# # plt.ylabel(r"$\epsilon$")

# # print(temperature)

# # Problem 9 

# ind_40 = np.argmax(cv40)
# ind_60 = np.argmax(cv60)
# ind_80 = np.argmax(cv80)
# ind_100 = np.argmax(cv100)

# tc_40 = temperature[ind_40]
# tc_60 = temperature[ind_60]
# tc_80 = temperature[ind_80]
# tc_100 = temperature[ind_100]

# tc_arr = np.array([tc_40, tc_60, tc_80, tc_100])
# L_arr = np.array([1/40, 1/60, 1/80, 1/100])

# res = linregress(L_arr, tc_arr)
# print(res.intercept, res.slope)
# print(res.intercept_stderr)
# plt.figure()
# plt.plot(L_arr, tc_arr, 'ko')
# plt.plot(L_arr, res.intercept + res.slope*L_arr, 'r', label='fitted line')


plt.show()


