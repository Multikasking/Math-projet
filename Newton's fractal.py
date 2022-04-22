# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 09:34:53 2022

@author: Raph
"""
import numpy as np
import random as rdn
import matplotlib.pyplot as plt

def P(z):
    return z**3-1
def newton_complexe(f,eps,h_range,w_range,itemax):
    y, x = np.ogrid[-1: 1: h_range*1j, -1: 1: w_range*1j]
    a_array=x+y*1j
    racine=np.zeros(a_array.shape)
    for h in range (h_range):
        for w in range (w_range):
                z=a_array[h][w]
                z=(2*z**3+1)/(3*z**2)
                for i in range (itemax):
                    z_0=z
                    z=(2*z**3+1)/(3*z**2)
                    if abs(z_0-z)<eps:
                        if np.imag(z)<=eps and np.imag(z)>=-eps:
                            racine[h][w]=30
                        elif np.imag(z)>eps:
                            racine[h][w]=60
                        elif np.imag(z)<eps:
                            racine[h][w]=90
                        break
               
    return racine
a=newton_complexe(P,10**-6,1800,1800,100)
plt.figure(figsize=(20,20))
plt.axis("off")
plt.imshow(a, cmap='plasma')

