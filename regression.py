# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Daniil Huryn)s
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sp
from scipy.optimize import curve_fit


def bacteria(num_bacteria_and_nutrient, g_max, K_const, a_const, Mu_const):  # the differential equations for the model

    dndt = [
        num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] / (num_bacteria_and_nutrient[1] + K_const) -
        num_bacteria_and_nutrient[0] * Mu_const,
        -a_const * num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] / (
                    num_bacteria_and_nutrient[1] + K_const)]  # differential equations

    return dndt  # return differential equation


df = pd.read_excel(r'D:\newDownload\Bacteria_data.xlsx', sheet_name='Sheet1')
data = df.values
t = data[:, 0]


##LogFit to find mu 
HighBactData = np.log(data[31:])

t2 = t[31:]
t2 = np.vstack([np.ones(t2.size), t2]).T
a2, mu_guess = (np.linalg.lstsq(t2, HighBactData[:, 1], rcond=None)[0])

logfitfor_mu = [a2 + mu_guess * i for i in data[32:, 0]]
exp_fit_for_mu = np.exp(logfitfor_mu)

###rofactor 
a_guess = np.amax(data) / 0.2
nutrient_amount = [0.2 - i / a_guess for i in data[0:52, 1]]
for i in np.arange(42,52):
    nutrient_amount[i] = 0

#####LogFit to find gmax
low_bact_data = np.log(data[0:29])

t1 = t[0:29]
t1 = np.vstack([np.ones(t1.size), t1]).T
a, g_max_guess = (np.linalg.lstsq(t1, low_bact_data[:, 1], rcond=None)[0])
g_max_guess = g_max_guess - mu_guess
logfit_of_g_max = [a + g_max_guess * i for i in data[0:29, 0]]
expfit_of_g_max = np.exp(logfit_of_g_max)

#### Find K                 TO-DO
yk = np.log(data[20:41])

rough_deriv_bact = [(data[i+1][1]- data[i][1]) for i in np.arange(0,51)]
rough_deriv_bact+=[rough_deriv_bact[50]]     #since we only have 51 differences in a 52 sized array
K_guess_list = [(data[i][1] * g_max_guess * nutrient_amount[i] / (data[i][1] * mu_guess + rough_deriv_bact[i]) - nutrient_amount[i]) for i in np.arange(32,41)]

K_guess = np.average(K_guess_list)


input = np.hstack([data[1], nutrient_amount])

#v, k = curve_fit(bacteria, [data[1], nutrient_amount], t, g_max_const, 1, a_guess, mu_guess)[0]



# plt.plot(A2[:,1], y2res, 'o', label='data')
# plt.plot(data1[:,0], data1[0:28,1], 'o', label='data')
# plt.plot(data[:,0], data[:,1], 'o', label='data')
# plt.semilogy(data[:,0], data[:,1], 'o', label='data')
