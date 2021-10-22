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


# In[5]:


def Markov_Chain(alpha=1, beta=1, delta=1, gamma=1):
    P = np.array([[1-alpha, 0, alpha*beta, 0],
                  [0, 1-beta, alpha*beta, 0],
                  [0, gamma, 1-gamma, delta], 
                  [0,0,0,1]])
    #fix state 
    state=np.array([[1.0, 0.0, 1.0, 0.0]])
    stateHist=state
    dfStateHist=pd.DataFrame(state)
    distr_hist = [[0,0,0,0]]
    
    for x in range(25): #steady state 
        state=np.dot(state,P) 
        stateHist=np.append(stateHist,state,axis=0)
        dfDistrHist = pd.DataFrame(stateHist)
        plot = plt.plot(dfDistrHist), plt.xlabel('Time'), plt.ylabel('Colonies')

    
    return plot, P, dfDistrHist

interact(Markov_Chain,alpha=(0.000,1,0.001), beta=(0.000,1,0.001), delta=(0.000,1,0.001), gamma=(0.000,1,0.001))


# In[ ]:





# In[ ]:




