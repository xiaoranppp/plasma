# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 18:33:23 2018

@author: Xiaoran Peng
"""
import numpy as np
import matplotlib.pyplot as plt
def process2(constant1,electron_t,constantup,constantdown,density_list,energy):
    k=constant1*(electron_t**constantup)*np.exp(constantdown/electron_t)
    r=k
    for reactor in density_list:
        r=r*reactor
    return k,r,energy
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
def pressure(N,T):
    P=50*N*T/(Ninitial*300)
    return P*(10**-3)
def electron_temperature_change(electron_density,power,volume,rate_list,energy_list):
    target=power/volume-sum((np.array(rate_list)*np.array(energy_list)))
    #print(sum((np.array(rate_list)*np.array(energy_list))))
    target=target/(1.5*electron_density*kb_ev)
    #print(target)
    return target

N0=1.6*(10**19)
V=0.0015
T=300
Tion=300
cl1=10**9
kb_ev=1
volume=1.5*(10**(3))
Te=0.5
ne=10**9
unit=0.02

Mcl2=1.66*(10**-27)*71
Mcl=1.66*(10**-27)*35.5
me=9.1*(10**-31)

cl_lj=2*(10**-8)
cl2_lj=2*(10**-8)

cl2=1.6*(10**19)

kb=1.38*(10**-23)
e=1.6*(10**-19)
cl=0
clv=0
cl11=0
cl21=0

Tev=0.0258

cl2value=[cl2]
clvvalue=[clv]
cl21value=[cl21]

Telist=[Te]

clvalue=[cl]
cl1value=[cl1]
cl11value=[cl11]

nevalue=[10**9]

#B=0.5
Ninitial=1.6*(10**15)+10**9
powerlist=[]
dtel=[]
total=[]


F=500*4.479*(10**(17))
#Twall=300
#Tin=300
#P0=0.05
#Tlist=[300]
#Tionlist=[300]
#plist=[0.05]

dlist=[]



for i in range(10000):
    if i<100:
        power=(0.1+i*50/100)/(1.6*(10**-19))
        #print(power)
    else:
        power=50/(1.6*(10**-19))
    Dion2=661/0.52
    Dion1=661/0.36
    Dneu=8776.7
    Dneuv=6269
    k1,r1,E1=process2(1.04*(10**-13),Te,-0.29,-8.84,[ne,cl2],4)
    k2,r2,E2=process2(5.12*(10**-14),Te,0.48,-12.34,[ne,cl2],11.5)
    k3,r3,E3=process2(2.14*(10**-13),Te,-0.07,-25.26,[ne,cl2],15.5)
    k4,r4,E4=process2(2.27*(10**-16),Te,1.92,-21.26,[ne,cl2],0)
    k5,r5,E5=process2(3.43*(10**-15),Te,-1.18,-3.98,[ne,cl2],0)
    k6,r6,E6=process2(2.94*(10**-16),Te,0.19,-18.79,[ne,cl2],11.9)
    k7,r7,E7=process2(1.22*(10**-16),Te,-0.99,-0.40,[ne,cl2],0.07)
    k8,r8,E8=process2(9*(10**-14),Te,-0.5,0,[ne,cl21],0)
    k9,r9,E9=process2(3.17*(10**-14),Te,0.53,-13.29,[ne,cl],12.96)
    k10,r10,E10=process2(9.02*(10**-15),Te,-0.92,-4.88,[ne,cl11],3.6)
    k11,r11,E11=process2(3.62*(10**-15),Te,0.72,-25.38,[ne,cl11],0)
    k12,r12,E12=process5((5*(10**-14)),[ne,cl11])
    k13,r13,E13=process5((5*(10**-14)),[cl1,cl11])
    k14,r14,E14=process5((5.4*(10**-16)),[cl2,cl1])
    k15,r15,E15=process6(4/2.405,(Dion1)*(1+Te/Tev),[cl1])
    k16,r16,E16=process6(4/2.405,(Dion1)*(1+Te/Tev),[cl11])
    k17,r17,E17=process6(4/2.405,Dneu,[2*cl])
    k18,r18,E18=process6(4/2.405,Dneuv,[clv])
    k19,r19,E19=process6(4/2.405,(Dion1)*(1+Te/Tev),[cl21])
    tdepend=(T/300)**0.5
    klist=np.array([k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12*tdepend,k13*tdepend,k14,k15,k16,k17,k18,k19])
    rlist=np.array([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12*tdepend,r13*tdepend,r14,r15,r16,r17,r18,r19])
    Elist=np.array([E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E19])
    cl2list=np.array([-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,0,1,1,1,1])
    clvlist=np.array([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,-1,0])
    cl21list=np.array([0,1,0,0,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,-1])
    cllist=np.array([2,0,1,0,1,0,0,2,-1,1,0,3,2,1,1,0,-2,0,0])
    cl1list=np.array([0,0,1,2,0,1,0,0,1,0,1,0,-1,-1,-1,0,0,0,0])
    cl11list=np.array([0,0,0,0,1,1,0,0,0,-1,-1,-1,-1,0,0,-1,0,0,0])
    nelist=cl21list+cl1list-cl11list
    
    
    cl2in=F/volume
    Nneural=cl+cl2+clv
    Ntotal=cl2+cl1+cl+clv+cl11+cl21
    Ninitial=1.6*(10**19)
    clout=-(cl/Nneural)*(1+(Ntotal-Ninitial)/Ninitial)*(F/volume)
    cl2out=-(cl2/Nneural)*(1+(Ntotal-Ninitial)/Ninitial)*(F/volume)
    clvout=-(clv/Nneural)*(1+(Ntotal-Ninitial)/Ninitial)*(F/volume)
   
    
    
    cl2=cl2+(sum(rlist*cl2list)+cl2in+cl2out)*unit*(10**-5)
    clv=clv+(sum(rlist*clvlist)+clvout)*unit*(10**-5)
    cl21=cl21+sum(rlist*cl21list)*unit*(10**-5)
    #print(ar1,sum(rlist*ar1list)*unit*(10**-5))
    cl=cl+(sum(rlist*cllist)+clout)*unit*(10**-5)
    cl1=cl1+(sum(rlist*cl1list))*unit*(10**-5)
    cl11=cl11+(sum(rlist*cl11list))*unit*(10**-5)
    
    ne=ne+sum(rlist*nelist)*unit*(10**-5)
    
    nevalue.append(ne)
    cl2value.append(cl2)
    clvvalue.append(clv)
    cl21value.append(cl21)
    clvalue.append(cl)
    cl1value.append(cl1)
    cl11value.append(cl11)
    dTe=electron_temperature_change(ne,power,volume,rlist,Elist)
    dlist.append(dTe)
    #print(cl1,Te)
#    if i>1170 :
#        dTe=0
    Te=dTe*unit*(10**-5)+Te
   # print(dTe)
    Telist.append(Te)
    #print(r1*E1,r2*E2,r3*E3,r6*E6,r7*E7,r9*E9,r10*E10)
plt.plot(Telist)
  
    
    
    
# 