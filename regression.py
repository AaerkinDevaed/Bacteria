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


def bacteria(num_bacteria_and_nutrient, t, g_max, K_const, a_const, Mu_const):  # the differential equations for the model, num_bact[0] is the number of bact, [1] is the number of nutrients
    
        dndt = [
        num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] / (num_bacteria_and_nutrient[1] + K_const) -
        num_bacteria_and_nutrient[0] * Mu_const,
        -a_const * num_bacteria_and_nutrient[0] * g_max * num_bacteria_and_nutrient[1] /
                     (num_bacteria_and_nutrient[1] + K_const)
                     ]  # differential equations, first is dn/dt, second is dro/dt

        return dndt  # return differential equation
    
def bacterialag(num_bacteria_and_nutrient, t, g_max, K_const, a_const, Mu_const):  # the differential equations for the model during lag time (only deaths)
    
    dndt = [ - num_bacteria_and_nutrient[0] * Mu_const,
    0]  # differential equations, first is dn/dt, second is dro/dt

    return dndt  # return differential equation

def bacteriaplot(t, g_max, K_const, a_const, Mu_const):  #odeint of the first equation for curvefit
    b = ((odeint(bacteria, init_val, t, args=(g_max, K_const, a_const, Mu_const)))).flatten()
    for i, x in enumerate(b):        #so we don't take log of 0 or neg numbers, we make all values <=0 equal to 0.0001 (they don't make sense physically either)
        if x<=0:
            b[i]=0.0001            
    return np.log(b)

df = pd.read_excel(r'D:\newDownload\Bacteria_data.xlsx', sheet_name='Sheet1')
data = df.values
t = data[:, 0]  #the time values of the given data


##LogFit to find mu
HighBactData = np.log(data[31:])          #data with high n and low rho

t2 = t[31:]
t2 = np.vstack([np.ones(t2.size), t2]).T
a2, mu_guess = (np.linalg.lstsq(t2, HighBactData[:, 1], rcond=None)[0])  #linear fit to find a mu estimate

###Finding Ro (needed for how our odeint works)
a_guess = 0.2 / np.amax(data[:,1])                       #estimate for a, assuming no deaths
nutrient_amount = [0.2 - i * a_guess for i in data[0:52, 1]] #again, this rho assumes bacteria alive are the only ones that used nutrients, so doesn't account for deaths
for i in np.arange(42,52):
    nutrient_amount[i] = 0.0000001

#####LogFit to find gmax
low_bact_data = np.log(data[0:29])              #high rho low n data

t1 = t[0:29]
t1 = np.vstack([np.ones(t1.size), t1]).T 
a, g_max_guess = (np.linalg.lstsq(t1, low_bact_data[:, 1], rcond=None)[0]) #linear fit to estimate g_max, not accounting for mu
g_max_guess = g_max_guess - mu_guess                                       #account for mu

#### Find K                
yk = np.log(data[20:41])                                                   #middle data to find K, when rho is small but still a factor

rough_deriv_bact = [(data[i+1][1]- data[i][1]) for i in np.arange(0,51)] #only correct for when time difference is 1
rough_deriv_bact+=[rough_deriv_bact[50]]     #since we only have 51 differences in a 52 sized array, but fortunately it is close to 1 on average for 20-41, so good enough for our guess
K_guess_list = [(data[i][1] * g_max_guess * nutrient_amount[i] / (data[i][1] * mu_guess + rough_deriv_bact[i]) - nutrient_amount[i]) for i in np.arange(32,41)]

K_guess = np.average(K_guess_list)          #veeeery rough estimate, see report

data_with_ro = np.vstack([data[:,1], nutrient_amount]).T               

init_val = [data[0][1], 0.2]

pass_data=(np.array([np.log(data_with_ro[:,0]),np.log(data_with_ro[:,1])]).T).flatten()           #data we pass to curve_fit


g_final, K_final, a_final, mu_final = curve_fit(bacteriaplot, t, pass_data, p0=(g_max_guess, K_guess, a_guess, mu_guess))[0]      #actual curve_fitting for  parameters, main result of project
print(g_final, K_final, a_final, mu_final)


res = ((odeint(bacteria, init_val, t, args=(g_final, K_final, a_final, mu_final))))   #odeint generated curve with our parameters from curve_fit
avg_error = np.average([1 - res[i,0]/data[i,1] for i in np.arange(0, 52)])            #find average error to see goodness of fit (alongside r^2)

plt.semilogy(data[:,0], data[:,1], 'o', label='data')       #given data
plt.semilogy(data[:,0], res[:,0], 'o', label='fit')        #plot of fitted model
plt.ylabel('CFU / ml')                                  #put axis titles on graph
plt.xlabel('Time (h)')
plt.legend(loc='best')
r1 = np.corrcoef(res[:,0], data[:,1])[0,1]**2                 #r^2 value
plt.title("Data and Fitted model of Bacteria growth, R^2 = %f" %(r1))
plt.figure()    
modellag = ((odeint(bacterialag, init_val, np.arange(0,4.5), args=(g_final, K_final, a_final, mu_final)))) #asked for 150 hour model - lag time 
modelgrow = ((odeint(bacteria, modellag[4], np.arange(4.50001,150), args=(g_final, K_final, a_final, mu_final)))) #asked for 150 hour model - growth time
model = np.vstack([modellag, modelgrow])
plt.semilogy(np.hstack([np.arange(0,4.5), np.arange(4.50001,150)]), model[:,0], 'o', label='ODEINT')
plt.ylabel('CFU / ml (Log)')                                  #put axis titles on graph
plt.xlabel('Time (h) ')
plt.legend(loc='best')
plt.title('Model for 150 hours using odeint')

# Fitting a 3rd order polynomial.

plt.figure()

data_to_3rd = np.vstack([np.ones(len(data[:,0])), data[:,0], data[:,0]**2, data[:,0]**3]).T
e1, e2, e3, e4  = np.linalg.lstsq(data_to_3rd, data[:,1], rcond=None)[0] #curvefit for 3rd order poly.
r2 = np.corrcoef([e1+e2*data[i,0]+e3*data[i,0]**2+e4*data[i,0]**3 for i in np.arange(0,52)], data[:,1])[0,1]**2       #r^2
plt.plot(data[:,0], e1+e2*data[:,0]+e3*data[:,0]**2+e4*data[:,0]**3, label='Fitted line')
plt.plot(data[:,0], data[:,1], 'o', label='data')
plt.ylabel('CFU / ml')                                  #put axis titles on graph
plt.xlabel('Time (h) ')
plt.legend(loc='best')
plt.title('Polynomial Fit with 3rd order and Data, R^2 = %f' %r2)


# Fitting a 10th order polynomial.

plt.figure()

data_to_10th = np.vstack([np.ones(len(data[:,0])), data[:,0], data[:,0]**2, data[:,0]**3, data[:,0]**4, data[:,0]**5, data[:,0]**6, data[:,0]**7, data[:,0]**8, data[:,0]**9, data[:,0]**10]).T 
p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 = np.linalg.lstsq(data_to_10th, data[:,1], rcond=None)[0] #curvefit for 10th order poly.
#looks bad and unreadable, unfortunately couldn't find other quick way to make a 10th order polynomial

r3 = np.corrcoef([p1+p2*data[i,0]+p3*data[i,0]**2+p4*data[i,0]**3 + p5*data[i,0]**4+p6*data[i,0]**5+p7*data[i,0]**6+p8*data[i,0]**7+p9*data[i,0]**8+p10*data[i,0]**9+p11*data[i,0]**10 for i in np.arange(0,52)], data[:,1])[0,1]**2    #r^2
plt.plot(data[:,0], (p1+p2*data[:,0]+p3*data[:,0]**2+p4*data[:,0]**3+p5*data[:,0]**4+p6*data[:,0]**5+p7*data[:,0]**6+p8*data[:,0]**7+p9*data[:,0]**8+p10*data[:,0]**9+p11*data[:,0]**10), label='Fitted line')
plt.plot(data[:,0], data[:,1], 'o', label='data')
plt.ylabel('CFU / ml')                                  #put axis titles on graph
plt.xlabel('Time (h) ')
plt.legend(loc='best')
plt.title('Polynomial Fit with 10th order and Data, R^2 = %f' %r3)