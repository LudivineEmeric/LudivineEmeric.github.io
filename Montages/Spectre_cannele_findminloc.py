#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:07:34 2021

code pour trouver les annulations dans la manip du spectre cannelé

@author: ludivine emeric
"""

import numpy as np
#from math import pi as pi
import matplotlib.pyplot as plt

lambdamin,lambdamax=400,800

#dossier="/home/ludivine.emeric/Documents/Python/MIMnanogap/SIMU/20180305_8021_conic0-90_12-24_gap21A_w248_n58/pluslbda/"
#  reflecTM=np.empty([len(tablambda), len(tabphi)],dtype=float)
#  	fname=dossier+"ReflecMoy_MIMnanogap99211_thetaphi_de"+str(tabtheta.min())+"deg_a_"+str(tabtheta.max())+"deg_pas_de_"+str(pas)+"_w"+str(wmim*1e3)+"_per"+str(per*1e3)+"_gap"+str(hALD*1e3)+"_nALD"+str(nALD*1e2)+".txt"
#18

dossier = "/Users/ludivine/Desktop/"
fname   = dossier+"spectre_cannele.txt"

text = open(fname, "r").readlines()
final=[]
i=0
for line in text:
    line_new=''
    for a,b in enumerate(line):
        if b==',': line_new=line_new+"."
        else: line_new = line_new+b
    final.append(line_new)
    i=i+1

        
c=len(text)
tab=np.loadtxt(final, skiprows=17, delimiter="\t", max_rows=c-17-1)
           
tablambda  = tab[:,0]
spectre = tab[:,1]
  
 
from scipy.signal import argrelextrema
index=argrelextrema(spectre,np.less)
index_interf=np.empty([0],dtype=float)
mask1=tablambda<lambdamax
mask2=tablambda>lambdamin
mask=mask1 & mask2
ll=tablambda[mask]
ss=spectre[mask]
tt=np.empty(len(ss),dtype=float)


plt.plot(ll,ss)

N=len(ll)

ss_min=np.min(ss)
# critere pour la limite :
inter_lim=2*ss_min

tt=ss
## lisser le spectre pour ameliorer reperage ?
for j in range(2,N-3):
#    tt[j]=(ss[j+2]+ss[j+1]+ss[j]+ss[j-1]+ss[j-2])/5
    tt[j]=(ss[j+1]+ss[j]+ss[j-1])/3
    
plt.plot(ll,ss,'r',ll,tt,'k')

derss=np.empty([N-1],dtype=float)
dertt=np.empty([N-1],dtype=float)



for j in range(0,N-2):
    dertt[j]=(tt[j+1]-tt[j])/(ll[j+1]-ll[j])
    
    
derss=dertt
## lisser la derive ??
# for j in range(2,N-3):
#     derss[j]=(dertt[j+2]+dertt[j+1]+dertt[j]+dertt[j-1]+dertt[j-2])/5
    
for j in range(1,N-1):
    a=derss[j-1]<0
    b=derss[j]>0
    if a & b:
        if ss[j]<inter_lim:
            c=[j]
            index_interf=np.append(index_interf,c)

#plt.plot(ll[0:N-1],derss,'b')
plt.xlim(400,800)
#plt.ylim(-20000,70000)


vector = np.vectorize(np.int)        
index_interf2=vector(index_interf)
#print(len(index_interf2))
X,Y=ll[index_interf2],ss[index_interf2]
plt.plot(X,Y,'or')
Nx=len(X)
W=np.arange(Nx,1,-1)

fname_new="spectre_cannele_annulations.txt"
np.savetxt(fname_new,np.c_[W, X, Y],fmt='%10.4f')
