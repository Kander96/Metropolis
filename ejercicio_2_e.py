#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:43:11 2019

@author: kevin
"""
import numpy as np
from matplotlib import pyplot as plt

def promediar_magnetizacion(datos):
    n=len(datos)
    m=np.zeros(n)
    for j in range(n):
        m[j]=np.mean(datos[j,2:])
    
    return m

B=[-9.00,-3.00,-1.00,0.00,1.00,3.00,9.00]
B=[-1,-0.5,0,0.5,1.0]

for i in range(len(B)):
    datos=np.loadtxt('ejercicio_2_e_B_%.2f.txt' %B[i])
    
    m=promediar_magnetizacion(datos)
    
    plt.plot(1/datos[:,1],m,'.',label='B=%.1f' %B[i])

plt.grid()
plt.legend()
plt.xlabel('1/J*')
plt.ylabel('Magnetización por sitio')
plt.savefig('primer_aproach_frustacion.pdf')
plt.show()

B=-1
datos=np.loadtxt('ejercicio_2_e_B_%.2f.txt' %B)
m=promediar_magnetizacion(datos)

plt.plot(1/datos[:,1],m,'.',label='B=%.1f' %B)
plt.grid()
plt.legend()
plt.xlabel('1/J*')
plt.ylabel('Magnetización por sitio')
#plt.savefig('antiferro_B_0.pdf')
plt.show()
