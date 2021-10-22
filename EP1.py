#!/usr/bin/env python
# coding: utf-8

# In[1]:


##Imports 
#slider 
#!pip install seir
#!pip install ipywidgets
#!jupyter nbextension enable --py widgetsnbextension
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from IPython.display import display, clear_output, set_matplotlib_formats
set_matplotlib_formats('pdf', 'svg')


# In[13]:


def end_point(Donors, Recipients, Transconjugants, C, e, psi, gamma):
# slider for input values
    d = Donors 
    r = Recipients
    t = Transconjugants 
    
    R_dot = psi*C*r - gamma*C*r*(d + t)
    D_dot = psi*C*d
    T_dot = psi*C*t + gamma*C*r*(d + t)
    C_dot = -psi*C*(r + d + t)* e 
    
    X = ['Recipients','Donors','Transconjugants','Resource Concentration']
    Y = [R_dot, D_dot, T_dot, C_dot]
  
    X_axis = np.arange(len(X))
  
    plt.bar(X_axis - 0.4, Y, 0.4, label = 'Cells / mL')

    plt.xticks(X_axis, X)
    plt.xlabel("Strains")
    plt.ylabel("Change in Cell Count")
    plt.title("Change in Cell Density")
    plt.legend()
    plt.show()
    
interact(end_point,
         Donors=(0.000,100,1), 
         Recipients=(0.000,100,1), 
         Transconjugants=(0.000,100,1), 
         C=(0.000,1,0.001),
         e=(0.000,1,0.001),
         psi=(0.000,1,0.001), 
         gamma=(0.000,1,0.001))


# In[3]:


def end_point_monod(psi, gamma, Q, C):
    psi_max = ((C)/(Q+C))/psi
    gamma_max = ((C)/(Q+C))/gamma
    
    X = ['Psi_max', 'Gamma_max']
    Y = [psi_max, gamma_max]
  
    X_axis = np.arange(len(X))
  
    plt.bar(X_axis - 0.2, Y, 0.4, label = 'Maximum Rates at Half Saturation Concentration')

    plt.xticks(X_axis, X)
    plt.xlabel("Plasmid Growth and Transfer Rates")
    plt.ylabel("Maximum Rate Threshold")
    plt.title("Monod Function of Resource Concentration")
    plt.legend()
    plt.show()
    
interact(end_point_monod,
         Q=(0.000,100,1), 
         C=(0.000,1,0.001),
         psi=(0.000,1,0.001), 
         gamma=(0.000,1,0.001))


# In[4]:


def end_point_time(Donors, Recipients, Transconjugants, C, psi, gamma):
    d = Donors 
    r = Recipients
    t = Transconjugants 
    
    N = r + d + t 
    Y = t / r 
    S = d / N
    beta = gamma / psi 
    
    R_dot = psi*C*r - gamma*C*r*(d + t)
    D_dot = psi*C*d
    T_dot = psi*C*t + gamma*C*r*(d + t)    
    N_dot = psi*N
    
    R_prime = R_dot / N_dot
    D_prime = D_dot / N_dot
    T_prime = T_dot / N_dot
    
    X = ['Recipients','Donors','Transconjugants']
    Y = [R_prime, D_prime, T_prime]
  
    X_axis = np.arange(len(X))
  
    plt.bar(X_axis - 0.2, Y, 0.4, label = 'Cells / mL')

    plt.xticks(X_axis, X)
    plt.xlabel("Strains")
    plt.ylabel("Change in Cell Count")
    plt.title("Rates of Transfomational Change with Respect to N")
    plt.legend()
    plt.show()

    return R_prime, D_prime, T_prime

interact(end_point_time,
         Donors=(0.000,100,1), 
         Recipients=(0.000,100,1), 
         Transconjugants=(0.000,100,1), 
         C=(0.000,1,0.001),
         psi=(0.000,1,0.001), 
         gamma=(0.000,1,0.001))


# In[5]:


from numpy import log as ln

def end_point_OD(Donors, Recipients, Transconjugants, N_0, ODa, ODb, delta_t):
    d = Donors 
    r = Recipients
    t = Transconjugants 
    N = r + d + t 
    
    psi = (np.log((ODb / ODa)))/delta_t
    gamma = psi*np.log(1+ (t /r)*(N /d)*(1/ (N -N_0)))
    
    X = ['Psi', 'Gamma']
    Y = [psi, gamma]
  
    X_axis = np.arange(len(X))
  
    plt.bar(X_axis - 0.2, Y, 0.4, label = 'Maximum Rates at Half Saturation Concentration')

    plt.xticks(X_axis, X)
    plt.xlabel("Plasmid Growth and Transfer Rates")
    plt.ylabel("Rate Threshold")
    plt.title("OD600 Function of Plasmid Transfer")
    plt.legend()
    plt.show()
    
    return psi, gamma


interact(end_point_OD,
         Donors=(0.000,100,1), 
         Recipients=(0.000,100,1), 
         Transconjugants=(0.000,100,1), 
         N_0=(0.000,100,1),
         N=(0.000,100,1), 
         ODa=(0.000,2,0.001),
         ODb=(0.000,2,0.001), 
         delta_t = (0.000,5,1))


# In[19]:


#ABM / Markov Addition 
def EP_Markov(Donors, Recipients, Transconjugants, N_0, ODa, ODb, delta_t):
    d = Donors 
    r = Recipients
    t = Transconjugants 
    N = r + d + t 
    
    psi = (np.log((ODb / ODa)))/delta_t #growth rate
    gamma = psi*np.log(1+ (t /r)*(N /d)*(1/ (N - N_0))) #transfer rate 
    
    P = np.array([[1-psi, 0, psi*gamma, 0],
                  [0, 1-psi, psi*gamma, 0],
                  [0, 1-gamma, 1-gamma, gamma], 
                  [0,0,0,1]])
    #fix state 
    state=np.array([[Donors, Recipients, 1.0, Transconjugants]]) #added one conjugated cell
    stateHist=state
    dfStateHist=pd.DataFrame(state)
    distr_hist = [[0,0,0,0]]
    
    for x in range(25): #steady state 
        state=np.dot(state,P) 
        stateHist=np.append(stateHist,state,axis=0)
        dfDistrHist = pd.DataFrame(stateHist)
        plot = plt.plot(dfDistrHist), plt.xlabel('Time'), plt.ylabel('Colonies')

    
    return plot, P, dfDistrHist
    #transconjugants --> blue
    #donors --> pink 
    #Recipients --> brown
    #Conjugated cells --> green 

interact(EP_Markov,
         Donors=(0.000,100,1), 
         Recipients=(0.000,100,1), 
         Transconjugants=(0.000,100,1), 
         N_0=(0.000,100,1),
         N=(0.000,100,1), 
         ODa=(0.000,2,0.001),
         ODb=(0.000,2,0.001), 
         delta_t = (0.000,5,1))


# In[ ]:




