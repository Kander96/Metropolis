# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:43:55 2019

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


paso=np.linspace(0.1,5,50)
a=[]

for l in paso:
    datos=np.loadtxt('distribucion_gaussiana_delta_%.2f.txt' %l)
    
    
    #plt.hist(datos,bins=50)
    #plt.show()
    
    #np.mean(datos)
    #np.var(datos)
         
    aceptado=1
    rechazado=0
    for i in range(len(datos)-1):
        if(datos[i+1]==datos[i]):
            rechazado=rechazado+1
        else:
            aceptado=aceptado+1
    
    a.append(aceptado/len(datos)*100)
    #print(r'Para $\delta=$',l,', el porcentaje de aceptacion es:',aceptado/len(datos)*100,'%\n')
    #print(r'Para $\delta=$',l,', el porcentaje de rechazo es:',rechazado/len(datos)*100,'%\n')
plt.plot(paso,a,'.')
plt.xlabel('Paso de exploración')
plt.ylabel('Porcentaje de aceptación')
plt.grid()
plt.savefig('Porcentaje_de_aceptacion.pdf')
plt.show()

paso=[0.5,1.0,3.0]
for l in paso:    
    datos=np.loadtxt('distribucion_gaussiana_delta_%.2f.txt' %l)
    L=1000
    N=100
    M=int(len(datos)/L)-1
    ro_k=np.zeros(N)
    for j in range(M):
        for i in range(N):
            ro_k[i]=ro_k[i]+correlacion(datos[j*L:],i,L)
    ro_k=ro_k/M
    
    k=np.linspace(0,N-1,N)
    plt.plot(k,ro_k,'.',label=r'$\delta=%.1f$' %l)
    #plt.title(r'$\delta=%.2f$' %paso)
plt.xlabel('Distancia de correlacion $k$')
plt.ylabel(r'Correlación $\rho(k)$')
plt.legend()    
plt.grid()    
plt.savefig('Correlación vs delta.pdf')
plt.show()
