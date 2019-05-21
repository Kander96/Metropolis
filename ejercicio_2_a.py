# -*- coding: utf-8 -*-
"""
Created on Tue May 14 17:35:09 2019

@author: kevin
"""
import numpy as np
from matplotlib import pyplot as plt

#%%
T=[0.5,50.0,100.0,500.0,1000.0,5000.0,10000.0]
M=np.zeros(len(T))
for i in range(len(T)):
    datos=np.loadtxt('ejercicio_2_a_T_%.1f.txt' %T[i])
    
    m=np.zeros(10000)
    for j in range(len(datos)):
        m=m+datos[j,:]
    m=m/len(datos)
    
    x=np.linspace(1,len(m),len(m))
    
    plt.plot(x,m,label='T=%.1f' %T[i])
    
    M[i]=np.mean(m[9000:])
plt.grid()
plt.legend()
plt.savefig('somethingsomething.pdf')
plt.show()

plt.plot(T,M,'o')
#plt.xscale('log')
plt.xlabel('Temperatura (K)')
plt.ylabel('magnetización por sitio')
plt.grid()
plt.xlim([-300,10300])
plt.ylim([-1.05,0.05])
plt.savefig('TvsM.pdf')
plt.show()
#%%

datos=np.loadtxt('ejercicio_2_a.txt')

x=np.linspace(min(datos[:,0]),max(datos[:,0]),100)
y=np.tanh(x)

N=len(datos)
M=np.zeros(N)
for i in range(N):
    M[i]=np.mean(datos[i,1:])

plt.plot(datos[:,0],M,'.',label='$<s_i>$')
plt.plot(x,y,label='$Tanh(B*)$')
plt.xlabel('B*')
plt.ylabel('Magnetización por sitio')
plt.ylim(-1.2,1.2)
plt.grid()
plt.legend(loc='upper left')
plt.savefig('TvsB_T.pdf')
plt.show()
