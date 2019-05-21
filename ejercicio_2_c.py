# -*- coding: utf-8 -*-
"""
Created on Sat May 18 20:10:22 2019

@author: kevin
"""
import numpy as np
from matplotlib import pyplot as plt

dim=[4,8,16,32,64]
for j in dim:
    datos=np.loadtxt('ejercicio_2_c_dim_%i.txt' %j)
    
    N=len(datos)
    M=np.zeros(N)
    for i in range(N):
        M[i]=np.mean(abs(datos[i,1:]))
    
    x=np.linspace(0,3,100)
    #y=(1-(np.sinh(2/x))**(-4))**(1/8)
    
    plt.plot(1/datos[:,0],M,'.',label='L=%i' %j)
    #plt.plot(x,y,'r')
plt.grid()
plt.legend()
plt.show()