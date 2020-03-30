# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Daniil Huryn)s
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint


def bacteria(num_bacteria_and_nutrient, t, g_max, K_const, a_const, Mu_const):  # the differential equations for the model
    
        dndt = [
        num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] / (num_bacteria_and_nutrient[1] + K_const) -
        num_bacteria_and_nutrient[0] * Mu_const,
        -a_const * num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] /
                     (num_bacteria_and_nutrient[1] + K_const)]  # differential equations

        return dndt  # return differential equation
    
def bacterialag(num_bacteria_and_nutrient, t, g_max, K_const, a_const, Mu_const):  # the differential equations for the model during lag time (only deaths)
    
    dndt = [ - num_bacteria_and_nutrient[0] * Mu_const,
    0]  # differential equations

    return dndt  # return differential equation

def bacteriaplot(t, g_max, K_const, a_const, Mu_const):
    b = ((odeint(bacteria, init_val, t, args=(g_max, K_const, a_const, Mu_const)))).flatten()
    for i, x in enumerate(b):
        if x<=0:
            b[i]=0.0001
    return np.log(b)

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
a_guess = 0.2 / np.amax(data)
nutrient_amount = [0.2 - i * a_guess for i in data[0:52, 1]]
for i in np.arange(42,52):
    nutrient_amount[i] = 0.0000001

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

data_with_ro = np.vstack([data[:,1], nutrient_amount]).T

init_val = [data[0][1], 0.2]

pass_data=(np.array([np.log(data_with_ro[:,0]),np.log(data_with_ro[:,1])]).T).flatten()


g_final, K_final, a_final, mu_final = curve_fit(bacteriaplot, t, pass_data, p0=(g_max_guess, K_guess, a_guess, mu_guess))[0]



res = ((odeint(bacteria, init_val, t, args=(g_final, K_final, a_final, mu_final))))


plt.semilogy(data[:,0], data[:,1], 'o', label='data')       #given data
plt.semilogy(data[:,0], res[:,0], 'o', label='fit')        #plot of fitted model
plt.xlabel('CFU / ml')                                  #put axis titles on graph
plt.ylabel('Time (h)')
plt.legend(loc='best')
plt.title('Data and Fitted model of Bacteria growth')
plt.figure()    
modellag = ((odeint(bacteria, init_val, np.arange(0,4.5), args=(g_final, K_final, a_final, mu_final)))) #asked for 150 hour model - lag time 
modelgrow = ((odeint(bacteria, modellag[4], np.arange(4.50001,150), args=(g_final, K_final, a_final, mu_final)))) #asked for 150 hour model - growth time
model = np.vstack([modellag, modelgrow])
plt.semilogy(np.hstack([np.arange(0,4.5), np.arange(4.50001,150)]), model[:,0], 'o', label='ODEINT')
plt.xlabel('CFU / ml')                                  #put axis titles on graph
plt.ylabel('Time (h) (Log)')
plt.legend(loc='best')
plt.title('Model for 150 hours using odeint')
# plt.plot(data[:,0], data[:,1], 'o', label='data')
# plt.semilogy(data[:,0], data[:,1], 'o', label='data')
