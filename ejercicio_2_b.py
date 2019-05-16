# -*- coding: utf-8 -*-
"""
Created on Tue May 14 22:50:38 2019

@author: kevin
"""
import numpy as np
from matplotlib import pyplot as plt

def correlacion(datos,k,L):
    ro=(algo(datos,k,L)-np.mean(datos)**2)/(np.mean(datos**2)-np.mean(datos)**2)    
    return ro
    
def algo(datos,k,L):
    suma=0
    for i in range(L-k):
        suma=suma+datos[i]*datos[i+k]
    return suma/(L-k)


datos=np.loadtxt('ejercicio_2_b.txt')

x=np.linspace(1,len(datos),len(datos))

plt.plot(x,datos)

datos=datos[5000:]

L=10000
k=np.arange(0,2001,100)
N=len(k)
ro=np.zeros(N)
M=int(len(datos)/L)-1
for j in range(M):
    for i in range(N):
        ro[i]=ro[i]+correlacion(datos[j*L:(j+1)*L],k[i],L)
ro=ro/M


plt.plot(k,ro,'.')
plt.grid()
plt.ylabel('Correlación')
plt.xlabel('Distancia de correlación')
plt.savefig('Correlacion_2_b.pdf')
plt.show()