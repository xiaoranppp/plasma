# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:47:03 2018

@author: Xiaoran Peng
"""
import numpy as np
#elastic process
#pressure mtorr
import matplotlib.pyplot as plt
def process0(constant1,electron_t,constantup,density_list,gas_temperature):
    k=constant1*(electron_t**constantup)
    r=k
    for reactor in density_list:
        r=r*reactor
    energy_loss=(2*5.46*(10**-4)/32)*(3/2)
    energy_loss=energy_loss*(electron_t-gas_temperature*kb/e)
    return k,r,energy_loss
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
def neu_diffusion_coefficient(constant1,gas_density,reduce_mass,lj_radius,T):
    return (10**2)*(constant1/gas_density)*((8.67*(10**-23)*T/reduce_mass)**0.5)/(3.14*lj_radius*lj_radius)
def ion_diffusion_coefficient(ion_mobility,gas_density_1atm_gas_density,T):
    return ion_mobility*(gas_density_1atm_gas_density)*T*kb/e
def electron_temperature_change(electron_density,power,volume,rate_list,energy_list):
    target=power/volume-sum((np.array(rate_list)*np.array(energy_list)))
   # print(sum((np.array(rate_list)*np.array(energy_list))))
    target=target/(1.5*electron_density*kb_ev)
    return target
def inputrate(f1,F):
    return f1*F/volume
def outputrate(g,Ntotal,F):
    return -g*(1+(Ntotal-N0)/N0)*(F/volume)
def pressure(N,T):
    P=50*N*T/(Ninitial*300)
    return P*(10**-3)
def ki(m1,m2,l1,l2,T):
    tmp=(m1+(m1/2+m2/2))*3.14*kb*T
    tp=(tmp/(2*m1*(m1/2+m2/2)))
    return (1.5*kb*(25/32)*((tp)**0.5))/(3.14*((l1/2+l2/2)**2))*100
#me=....
N0=1.6*(10**15)
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
Mo2=1.66*(10**-27)*32
Mo=1.66*(10**-27)*16
me=9.1*(10**-31)

Ar_lj=3.542*(10**-8)
o2_lj=3.47*(10**-8)
o_lj=3.05*(10**-8)

Ar=1.6*(10**15)*0.8
o2=1.6*(10**15)*0.2
kb=1.38*(10**-23)
e=1.6*(10**-19)
a_r=0
o2v=0
o_2=0
o21=0
o=0
o1=0
Tev=0.0258
o11=0#o+
#2_pi_kb=8.66*(10**-23)
Arvalue=[Ar]
ar1value=[ar1]
a_rvalue=[a_r]
Telist=[Te]
o2value=[o2]
o2vvalue=[0]
o_2value=[0]
o21value=[0]
ovalue=[0]
o1value=[0]
o11value=[0]
nevalue=[10**9]
B=1
Ninitial=1.6*(10**15)+10**9
powerlist=[]
dtel=[]
total=[]
kk=[]

part1list=[]
part2list=[]
part3list=[]
part4list=[]
part5list=[]

F=500*4.479*(10**(17))
Twall=300
Tin=300
P0=0.05
Tlist=[300]
Tionlist=[300]
plist=[0.05]
for i in range(5000):
#    if i<100:
#        power=(0.1+i*500/100)/(1.6*(10**-19))
#    elif i>=100 and i<=105:
#        power=(500.0001+500.0001*(100-i)/5)/(1.6*(10**-19))
#    elif i>105 and i<2500:
#        power=0.0001
        
        #print(power)
    
#    
#    step=i%2500
#    if step<=100:
#          power=(0.0001+step*(1333.7001-0.0001)/100)/(1.6*(10**-19))
#    elif step>100 and step<=475:
#          power=1333.7001/(1.6*(10**-19))
#    elif step>475 and step<=575:
#          power=(1333.7001+1333.7001*(475-step)/100)/(1.6*(10**-19))
#    else:
#        power=0.0001/(1.6*(10**-19))
    
    if i<100:
        power=(0.1+i*200/100)/(1.6*(10**-19))
        #print(power)
    else:
        power=200/(1.6*(10**-19))
#    
    powerlist.append(power)
    Ntotal=Ar+a_r+ar1+o2+o2v+o_2+o21+o+o1+o11
    total.append(Ntotal)
    #Dneo=1850
    #
    
    #P=Ntotal/Ninitial*(50*10**-3)
    P=pressure(Ntotal,T)
    Dneo=neu_diffusion_coefficient(3/16,Ntotal,Mar/2,Ar_lj,T)#neu_diffusion_coefficient(3/16,Ntotal,Mar/2,Ar_lj)#2639
    Dneoo_2=neu_diffusion_coefficient(3/16,Ntotal,Mo2/2,o2_lj,T)#2906,neu_diffusion_coefficient(3/16,Ntotal,Mo2/2,o2_lj)    Dion=47000/21#ion_diffusion_coefficient(1.53,760/P)

    Dneoo2v=neu_diffusion_coefficient(3/16,Ntotal,Mo2/2,o2_lj,T)#2906neu_diffusion_coefficient(3/16,Ntotal,Mo2/2,o2_lj)
    Dneoo=neu_diffusion_coefficient(3/16,Ntotal,Mo/2,o_lj,T)#4683neu_diffusion_coefficient(3/16,Ntotal,Mo/2,o_lj)
    Dion=ion_diffusion_coefficient(1.53,760*T/(P*273),Tion)#659
    Dion2=ion_diffusion_coefficient(2.57,760*T/(P*273),Tion)#1106ion_diffusion_coefficient(2.57,760/P)
    Dion3=ion_diffusion_coefficient(2.57,760*T/(P*273),Tion)
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
    k12,r12,E12=process0((4.79*(10**-8)),Te,0.5,[ne,o2],T)
    k13,r13,E13=process2((6.86*(10**-9)),Te,0,-6.29,[ne,o2],6.29)
    k14,r14,E14=process2((3.49*(10**-9)),Te,0,-5.92,[ne,o2],5.92)
    k15,r15,E15=process2((1.07*(10**-9)),Te,-1.39,-6.26,[ne,o2],0)
   # k15,r15,E15=process2((3.4*(10**-9)),Te,-1.6,-0.13,[ne,o2],0)
    k16,r16,E16=process2((2.8*(10**-9)),Te,0,-3.72,[ne,o2],0.19)
    k17,r17,E17=process2((1.28*(10**-9)),Te,0,-3.67,[ne,o2],0.385)
    k18,r18,E18=process2((1.37*(10**-9)),Te,0,-2.14,[ne,o2],0.977)
    k19,r19,E19=process2((2.34*(10**-9)),Te,1.03,-12.3,[ne,o2],12.3)
    k20,r20,E20=process2((2.8*(10**-9)),Te,0,-3.53,[ne,o2v],-0.19)
    k21,r21,E21=process2((2.34*(10**-9)),Te,1.03,-12.11,[ne,o2v],12.11)
    k22,r22,E22=process2((6.86*(10**-9)),Te,0,-6.1,[ne,o2v],6.1)
    k23,r23,E23=process2((3.49*(10**-9)),Te,0,-5.73,[ne,o2v],5.73)
    k24,r24,E24=process2((1.07*(10**-9)),Te,-1.39,-6.26,[ne,o2v],0)
    #k24,r24,E24=process2((3.4*(10**-9)),Te,-1.6,-0.13,[ne,o2v],0)
    k25,r25,E25=process2((2.06*(10**-9)),Te,0,-1.163,[ne,o_2],-0.977)
    k26,r26,E26=process2((2.34*(10**-9)),Te,1.03,-11.32,[ne,o_2],11.32)
    k27,r27,E27=process2((6.86*(10**-9)),Te,0,-5.31,[ne,o_2],5.31)
    k28,r28,E28=process2((3.49*(10**-9)),Te,0,-4.94,[ne,o_2],4.94)
    k29,r29,E29=process2((1.07*(10**-9)),Te,-1.39,-6.26,[ne,o_2],0)
    k30,r30,E30=process5((5*(10**-8)),[ar1,o1])
    k31,r31,E31=process5((5*(10**-8)),[o21,o1])
    k32,r32,E32=process4(2.2*(10**-8),Te,-0.5,[o21,ne])
    k33,r33,E33=process5((4.9*(10**-11)),[o2,ar1])
    k34,r34,E34=process5((1.0*(10**-10)),[o2,a_r])
    k35,r35,E35=process6(4/2.405,(Dion2)*(1+Te/Tev),[o21])#problem can be
    k36,r36,E36=process6(4/2.405,Dneoo_2,[o_2])
    k37,r37,E37=process6(4/2.405,Dneoo2v,[o2v])
    k38,r38,E38=process6(4/2.405,Dneoo,[o])#**2
    k39,r39,E39=process6(4/2.405,Dneoo,[o])
    
    k40,r40,E40=process5(1*(10**-7),[o,ne])
    k41,r41,E41=process2(9*(10**-9),Te,0.7,-13.62,[ne,o],13.62)
    k42,r42,E42=process5(4*(10**-8),[o1,o11])
    k43,r43,E43=process5(2.3*(10**-10),[o1,o])
    k44,r44,E44=process5(6.4*(10**-12),[o,ar1])
    k45,r45,E45=process5(4.1*(10**-11),[o,a_r])
    k46,r46,E46=process5(2*(10**-11),[o11,o2])
    k47,r47,E47=process2(5.47*(10**-8),Te,0.324,-2.98,[ne,o1],0)
    k48,r48,E48=process5(2*(10**-12),[o2v,Ar+a_r+o2+o_2+o])#??
    k49,r49,E49=process5(5.66*(10**-10),[ar1,Ar])
    k50,r50,E50=process5(1*(10**-9),[o21,o2])
    k51,r51,E51=process5(1*(10**-9),[o11,o])
    k52,r52,E52=process6(4/2.405,(Dion3)*(1+Te/Tev),[o11])
    
    tdepend=(T/300)**0.5
    klist=np.array([k1,k2,k3,k4,k5,k6,k7,k8,k9*tdepend,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30*((300/T)**0.5),k31*((300/T)**0.5),k32,k33*((300/T)**0.78),k34*tdepend,k35,k36,k37,k38,k39,k40,k41,k42*((300/T)**-0.43),k43*tdepend,k44*tdepend,k45*tdepend,k46*tdepend,k47,k48*tdepend,k49*tdepend,k50*tdepend,k51*tdepend,k52])
    rlist=np.array([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r40,r41,r42,r43,r44,r45,r46,r47,r48,r49,r50,r51,r52])
    Elist=np.array([E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,E21,E22,E23,E24,E25,E26,E27,E28,E29,E30,E31,E32,E33,E34,E35,E36,E37,E38,E39,E40,E41,E42,E43,E44,E45,E46,E47,E48,E49,E50,E51,E52])
    
#    print(klist,i)
    #print(rlist,i)
 #   print(Elist,i)
    arreaction=[0,0,0,0,0,0,0,0,0,0,0]
    #o2
    o2reaction=[0,-1,-1,-1,-1,-1,-1,-1,1,0,0,0,0,1,0,0,0,0,0,1,0,-1,-1,1,1,1,B/2,0]
    #o
    oreaction=[0,2,2,1,0,0,0,0,0,0,2,2,1,0,0,2,2,1,1,1,2,0,2,0,0,0,-B,0]
    #o2(v)
    o2vreaction=[0,0,0,0,1,1,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0]
    #o-
    o1reaction=[0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0]
    #o2*
    o_2reaction=[0,0,0,0,0,0,1,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,0,0,0]
    #o2+
    o21reaction=[0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,-1,-1,1,0,-1,0,0,0,0]
    noreaction=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    noreaction1=[0]*39
    
    Arlist=np.array([0,-1,-1,1,0,1,0,0,1,1,1]+noreaction+[1,0,0,1,1,0,0,0,0,0]+[0,0,0,0,1,1,0,0,0,0,0,0,0])
    a_rlist=np.array([0,1,0,-1,-1,-1,1,1,-2,0,-1]+noreaction+[0,0,0,0,-1,0,0,0,0,0]+[0,0,0,0,0,-1,0,0,0,0,0,0,0])
    ar1list=np.array([0,0,1,0,1,0,-1,-1,1,-1,0]+noreaction+[-1,0,0,-1,0,0,0,0,0,0]+[0,0,0,0,-1,0,0,0,0,0,0,0,0])
    

    o2list=np.array(arreaction+o2reaction+[0,0,0,1,0,0,-1,0,1,0,0,0,1])
    olist=np.array(arreaction+oreaction+[0,-1,2,-1,-1,0,1,1,0,0,0,0,0])
    o2vlist=np.array(arreaction+o2vreaction+[0,0,0,0,0,0,0,0,-1,0,0,0,0])
    o1list=np.array(arreaction+o1reaction+[0,0,-1,-1,0,0,0,-1,0,0,0,0,0])
    o_2list=np.array(arreaction+o_2reaction+[0]*13)
    o21list=np.array(arreaction+o21reaction+[0,0,0,0,0,0,1,0,0,0,0,0,0])
    o11list=np.array(noreaction1+[0,1,-1,0,1,0,-1,0,0,0,0,0,-2])
    nelist=ar1list+o21list+o11list-o1list
    
    #inflow and outflow
    Arin=0.9*F/volume
    o2in=0.1*F/volume
    Nneural=Ar+a_r+o2+o2v+o_2+o
    
    Ninitial0=1.6*(10**15)
    Arout=-(Ar/Nneural)*(1+(P-P0)/P0)*(F/volume)
    o2out=-(o2/Nneural)*(1+(P-P0)/P0)*(F/volume)
    oout=-(o/Nneural)*(1+(P-P0)/P0)*(F/volume)
    o2vout=-(o2v/Nneural)*(1+(P-P0)/P0)*(F/volume)
    o_2out=-(o_2/Nneural)*(1+(P-P0)/P0)*(F/volume)
    a_rout=-(a_r/Nneural)*(1+(P-P0)/P0)*(F/volume)
    kk.append(Arin+Arout)
    Ar=Ar+(sum(rlist*Arlist)+Arin+Arout)*unit*(10**-5)
    a_r=a_r+(sum(rlist*a_rlist)+a_rout)*unit*(10**-5)
    ar1=ar1+sum(rlist*ar1list)*unit*(10**-5)
    #print(ar1,sum(rlist*ar1list)*unit*(10**-5))
    o2=o2+(sum(rlist*o2list)+o2in+o2out)*unit*(10**-5)
    o=o+(sum(rlist*olist)+oout)*unit*(10**-5)
    o2v=o2v+(sum(rlist*o2vlist)+o2vout)*unit*(10**-5)
    o1=o1+sum(rlist*o1list)*unit*(10**-5)
    o_2=o_2+(sum(rlist*o_2list)+o_2out)*unit*(10**-5)
    o21=o21+sum(rlist*o21list)*unit*(10**-5)
    o11=o11+sum(rlist*o11list)*unit*(10**-5)
    #print(Te,power)
    #ne=ar1+o21-o1
    ne=ne+sum(rlist*nelist)*unit*(10**-5)
    if ne<10**6:
        ne=10**6
    nevalue.append(ne)
    Arvalue.append(Ar)
    a_rvalue.append(a_r)
    ar1value.append(ar1)
    o2value.append(o2)
    ovalue.append(o)
    o2vvalue.append(o2v)
    o1value.append(o1)
    o_2value.append(o_2)
    o21value.append(o21)
    o11value.append(o11)
    dTe=electron_temperature_change(ne,power,volume,rlist,Elist)
#    if i>1170 :
#        dTe=0
    Te=dTe*unit*(10**-5)+Te
    if Te<0.1:
        Te=0.1
    
    
#    if i>=3000:
#        Te=2.0
    Telist.append(Te)
    dtel.append(dTe)
    print(dTe)
    #print(rlist[len(rlist)-5:len(rlist)])
    
    #print(Te)
   # print(ne,ar1,o21,o1)
   #start T!!!
    sumr=sum([r30,r31,r33,r42,r44,r46,r49,r50,r51])
    part1=1.5*sumr*kb*(Tion-T)  
#    sumen=(r15+r24+r29)*-2.51#+(r13+r22+r27)*-1.03+(r14+r23+r28)*-0.82+r48*-0.3+r32*-7.1+r34*-6.5
#    part2=1.5*sumen*1.6*(10**-19)
    #a_r
    sumk=(1/Nneural)*(a_r/ki(Mar,Mar,Ar_lj,Ar_lj,T)+o2/ki(Mo2,Mar,o2_lj,Ar_lj,T)+o_2/ki(Mo2,Mar,o2_lj,Ar_lj,T)+o2v/ki(Mo2,Mar,o2_lj,Ar_lj,T)+o/ki(Mo,Mar,o_lj,Ar_lj,T)+Ar/ki(Mar,Mar,Ar_lj,Ar_lj,T))
    part3=(1/sumk)*(T-Twall)/(4/2.405)**2
    part4=(F/volume)*(kb*1.5)*(Tin-T*(1+(P-P0)/P0))
    part5=1.5*ne*2*me*(Te*e/kb-T)*(Ar*k1/Mar+o2*k12/Mo2+o*k40/Mo)*kb
    #print(part1,part2,part3,part4,part5)
    Tek=Te*e/kb
    Ea=(kb*Tek/(e*4/2.405))
    Nion=ne
    s = 760/273/(P/T)
    o21part=o21*Mo2*((2.57*s)**2)
    o11part=o11*Mo*((3*s)**2)
    o1part=o1*Mo*((3*s)**2)
    ar1part=ar1*Mar*((1.53*s)**2)
    
    
    Tion=T+(10**-4)*(1/3/kb)*(1/Nion)*(Ea**2)*sum([o21part,o11part,o1part,ar1part])
    T=T+((part1-part3+part4+part5)/(1.5*kb*Ntotal))*unit*(10**-5)
    part1list.append(part1)
    #part2list.append(part2)
    part3list.append(part3)
    part4list.append(part4)
    part5list.append(part5)
    Tlist.append(T)
    Tionlist.append(Tion)
    plist.append(P)
    
plt.plot(Telist)
#    
#    
#    
