#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 00:10:40 2021

@author: Shashaank
"""
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags
import bitfuncs
import closure
import enumdags
import math
from bitfuncs import val2nice
from bitfuncs import nice2val
from bitfuncs import subsetsof

def countBits(number):
     
     #log function in base 2 
     #take only integer part
     return int((math.log(number) /
                math.log(2)) + 1);

def arr(A):
    count=countBits(A)    
    a=[]
    while count>=0:
        if(A & (1<<count)):
           a.append(2**count)
        count=count-1  
    #print(a) 
    return a  

def parents(A,par):
    count=countBits(A)    
    parent=[]
    while count>=0:
        if(A & (1<<count)):
           parent.append(par[count])
        count=count-1  
    #print(parent) 
    return parent  
#W=0b101100
#par=[0b0,0b0,0b11,0b0,0b1110,0b10101]
#parents(W,par)    

par=[0b0,0b0,0b11,0b0,0b1110,0b10101]
obs=0b111100
X=0b100000
Y=0b1000
Z=0b110011
W=0b111010  
#Z=0b100
#W=0b10000  
  
def nondescendant(Z,W,obs,par):
    z,w,parent,n,summ=arr(Z),arr(W),parents(W,par),[],0
    strat=enumdags.dag_indeps(par)
    trips=closure.closure(strat)
    #print(z)
    #print(w)
    #print(parent)
    for i in range(len(z)):
        m=0
        for j in range(len(w)):
            if(z[i] in bitfuncs.subsetsof(parent[j])):
                m=m+0
            if(z[i] not in bitfuncs.subsetsof(parent[j])):
                 if(enumdags.implies(trips, closure.Triplet(z[i],w[j],parent[j]))==True):
                    m=m+0
                 if(enumdags.implies(trips, closure.Triplet(z[i],w[j],parent[j]))==False):
                    m=m+1
                    break
        n.append(m)        
    print(n)            
    for x in n:
        summ=summ+x               
    if(summ==0):
       print("No node in Z is descendant from W")
       return True
    else:
       print("Some node in Z is descendant from W") 
       return False            
            
nondescendant(Z, W, obs, par)    