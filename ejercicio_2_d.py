# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from matplotlib import pyplot as plt

def termalizacion(datos):
    N=len(datos)
    term=np.zeros([N,3])
    for i in range(N):
        term[i,0]=datos[i,0]
        term[i,1]=datos[i,1]
        term[i,2]=funcion(datos[i,2:])
    return term

def funcion(vector):
    m=np.mean(vector[int(len(vector)/2):])
    i=0
    while(abs(vector[i]-m)>0.05):
        i=i+1
    return i

def correlacion(datos,k,L):
    ro=(algo(datos,k,L)-np.mean(datos)**2)/(np.mean(datos**2)-np.mean(datos)**2)    
    return ro
    
def algo(datos,k,L):
    suma=0
    for i in range(L-k):
        suma=suma+datos[i]*datos[i+k]
    return suma/(L-k)

def promediar_magnetizacion(datos):
    n=len(datos)
    m=np.zeros(n)
    for j in range(n):
        m[j]=np.mean(datos[j,2:])
    
    return m


datos=np.loadtxt('ejercicio_2_d_termalizacion.txt')

x=np.linspace(1,len(datos[0,2:]),len(datos[0,2:]))

for i in range(len(datos)):
    plt.clf() 
    plt.plot(x,datos[i,2:])
    plt.title('J*=%.2f B*=%.2f' %(datos[i,1],datos[i,0]))
    plt.grid()
    plt.ylim([-1,1])
    plt.pause(0.5)

term=termalizacion(datos)

plt.plot(term[:,0],term[:,2],'.')
plt.grid()
plt.savefig('termalizacion_antiferro.pdf')
plt.show()

N=len(datos)
k=2000
ro=np.zeros([N,k+3])
for j in range(N):
    ro[j,0]=datos[j,0]
    ro[j,1]=datos[j,1]

B=[-9.00,-3.00,-1.00,0.00,1.00,3.00,9.00]

for i in range(len(B)):
    datos=np.loadtxt('ejercicio_2_d_B_%.2f.txt' %B[i])
    
    m=promediar_magnetizacion(datos)
    
    plt.plot(1/datos[:,1],m,'.',label='B=%.1f' %B[i])

plt.grid()
plt.legend()
plt.xlabel('1/J*')
plt.ylabel('Magnetización por sitio')
#plt.savefig('primer_aproach_antiferro.pdf')
plt.show()

B=0
datos=np.loadtxt('ejercicio_2_d_B_%.2f.txt' %B)
m=promediar_magnetizacion(datos)

plt.plot(1/datos[:,1],m,'.',label='B=%.1f' %B)
plt.grid()
plt.legend()
plt.xlabel('1/J*')
plt.ylabel('Magnetización por sitio')
plt.savefig('antiferro_B_0.pdf')
plt.show()

