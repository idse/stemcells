# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 10:11:18 2016

@author: idse
"""

import numpy as np
import matplotlib.pyplot as plt

num_sims = 5
N        = 1000

y_init = 0
t_init = 3
t_end  = 7

c_theta = 0.7
c_mu    = 1.5
c_sigma = 0.06

def mu(y, t): 
    return c_theta * (c_mu - y)
        
def sigma(y, t): 
    return c_sigma

dt   = float(t_end - t_init) / N
dW   = lambda dt: np.random.normal(loc = 0.0, scale = np.sqrt(dt))

t    = np.arange(t_init, t_end, dt)
y    = np.zeros(N)
y[0] = y_init

for i_sim in range(num_sims):
    for i in range(1, t.size):
        a = mu(y[i-1], (i-1) * dt)
        b = sigma(y[i-1], (i-1) * dt)
        y[i] = y[i-1] + a * dt + b * dW(dt)
    plt.plot(t, y)

plt.show()