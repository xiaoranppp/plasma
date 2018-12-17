# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 00:16:46 2018

@author: Xiaoran Peng
""
# -*- coding: utf-8 -*-
"""

import numpy as np
import matplotlib.pyplot as plt

def process1(constant1,electron_t,constantup,constantdown,density_list,mass,gas_temperature):
    k=constant1*(np.exp(constantup/(electron_t+constantdown)))
    r=k
    for reactor in density_list:
        r=r*reactor
    energy_loss=(2*5.46*(10**-4)/40)*(3/2)
    energy_loss=energy_loss*(electron_t-gas_temperature*kb/e)
    return k,r,energy_loss
#inelastic and super_lastic process
def process2(constant1,electron_t,constantup,constantdown,density_list,energy):
    k=constant1*(electron_t**constantup)*np.exp(constantdown/electron_t)
    r=k
    for reactor in density_list:
        r=r*reactor
    return k,r,energy
#mixed state radiation trap S-1
def process3(constant,Ar,a_r,energy=0):
    return constant*(Ar_0/Ar),constant*(Ar_0/Ar)*a_r,energy
#recombination
def process4(constant,electron_t,constantup,density_list,energy=0):
    k=constant*(electron_t**constantup)
    r=k
    for reactor in density_list:
        r=r*reactor
    return k,r,energy
#penning ionization
def process5(constant,density_list,energy=0):
    k=constant
    r=k
    for reactor in density_list:
        r=r*reactor
    return k,r,energy
# wall reaction
def process6(diffusion_len,diffusion_coef,density_list,energy=0):
    k=diffusion_coef/(diffusion_len**2)
    r=k
    for reactor in density_list:
        r=r*reactor
    return k,r,energy

def electron_temperature_change(electron_density,power,volume,rate_list,energy_list):
    target=power/volume-sum((np.array(rate_list)*np.array(energy_list)))
   # print(sum((np.array(rate_list)*np.array(energy_list))))
    target=target/(1.5*electron_density*kb_ev)
    return target

#me=....
Ninitial=1.6*(10**15)
Ar_0=(10**14)*3
V=0.0015
T=300
Tion=300
ar1=10**9
kb_ev=1
volume=1.5*(10**(3))
Te=0.5
ne=10**9
unit=0.02

Mar=1.66*(10**-27)*40
me=9.1*(10**-31)

Ar_lj=3.542*(10**-8)

Ar=1.6*(10**15)*0.8
kb=1.38*(10**-23)
e=1.6*(10**-19)
a_r=0

Tev=0.0258

#2_pi_kb=8.66*(10**-23)
Arvalue=[Ar]
ar1value=[ar1]
a_rvalue=[a_r]
Telist=[Te]
nevalue=[10**9]

dtel=[]
total=[]
kk=[]

F=0*4.479*(10**(17))
P0=0.05
Tlist=[300]

for i in range(5):
    if i<100:
        power=(0.1+i*200/100)/(1.6*(10**-19))
        #print(power)
    else:
        power=200/(1.6*(10**-19))
#    
    
    Ntotal=Ar+a_r+ar1
    total.append(Ntotal)
    
    Dneo=2633
    Dion=659
    k1,r1,E1=process1((3.9*(10**-7)),Te,-4.6,0.5,[ne,Ar],Mar,T)
    k2,r2,E2=process2(2.5*(10**-9),Te,0.74,-11.6,[ne,Ar],11.6)
    k3,r3,E3=process2(2.3*(10**-8),Te,0.68,-16,[ne,Ar],16)
    k4,r4,E4=process2(4.3*(10**-10),Te,0.74,0,[ne,a_r],-11.6)
    k5,r5,E5=process2(6.8*(10**-9),Te,0.67,-4.4,[ne,a_r],4.4)
    k6,r6,E6=process3(10**5,Ar,a_r)
    k7,r7,E7=process4(4.3*(10**-13),Te,-0.63,[ne,ar1])
    k8,r8,E8=process4(1.95*(10**-27),Te,-4.5,[ne,ne,ar1])
    k9,r9,E9=process5(1.2*(10**-9),[a_r,a_r])
    k10,r10,E10=process6(4/2.405,(Dion)*(1+Te/Tev),[ar1])
    k11,r11,E11=process6(4/2.405,Dneo,[a_r])
    #o
    
    
    tdepend=(T/300)**0.5
    klist=np.array([k1,k2,k3,k4,k5,k6,k7,k8,k9*tdepend,k10,k11])
    rlist=np.array([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11])
    Elist=np.array([E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11])


    Arlist=np.array([0,-1,-1,1,0,1,0,0,1,1,1])
    a_rlist=np.array([0,1,0,-1,-1,-1,1,1,-2,0,-1])
    ar1list=np.array([0,0,1,0,1,0,-1,-1,1,-1,0])

    nelist=ar1list
    
    #inflow and outflow
    Arin=F/volume
    Nneural=Ar+a_r
    
    Ntotal=Ar+a_r+ar1
    Arout=-(Ar/Nneural)*(1+(Ntotal-Ninitial)/Ninitial)*(F/volume)
    a_rout=-(a_r/Nneural)*(1+(Ntotal-Ninitial)/Ninitial)*(F/volume)
    Ar=Ar+(sum(rlist*Arlist)+Arin+Arout)*unit*(10**-5)
    a_r=a_r+(sum(rlist*a_rlist)+a_rout)*unit*(10**-5)
    ar1=ar1+sum(rlist*ar1list)*unit*(10**-5)
    #print(Te,power)
    #ne=ar1+o21-o1
    ne=ne+sum(rlist*nelist)*unit*(10**-5)
    if ne<10**6:
        ne=10**6
    nevalue.append(ne)
    Arvalue.append(Ar)
    a_rvalue.append(a_r)
    ar1value.append(ar1)
    dTe=electron_temperature_change(ne,power,volume,rlist,Elist)
#    if i>1170 :
#        dTe=0
    Te=dTe*unit*(10**-5)+Te
    Telist.append(Te)
    dtel.append(dTe)
    

plt.plot(Telist)
 
