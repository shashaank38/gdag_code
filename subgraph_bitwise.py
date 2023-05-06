#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 20:26:57 2021

@author: Shashaank
"""
alph=['A','B','C','D','E','F','G','H','I','J']
import bitfuncs
import closure
import enumdags
from bitfuncs import val2nice
from closure import Triplet
par=[0b0,0b1,0b11,0b0,0b0,0b11000, 0b100000,0b1010100]#array of parents
obs=0b11100110 #observed nodes
W="EFGB" #deleted nodes......this way it works even if you dont enter the nodes in alphabetical order
def subgraph_string(par,obs,W):
    w,subgraph,subgraph1=[],par,[]
    for i in W:
       w.append(bitfuncs.nice2val(i))
    print (w) 
    for i in range(len(w)):
       position= ord(W[i])-65
       subgraph[position]=0b0 # deleting parents llinks to the deleted node
       j=position
       while j<len(par):
           #if(w[i] in bitfuncs.subsetsof(par[j])):
               subgraph[j]=subgraph[j] & ~(1 << position)  #deleteing links from deleted node
               j=j+1
     #return subgraph

#program ends
       
     # I could change par directly in the above steps...just not doing incasse any issue arises later in other functions due to mutiple pars               
##To check quickly.. delete these lines later    
    print(subgraph) 
    
    for i in range(len(subgraph)):       
        subgraph1.append(bitfuncs.val2nice(subgraph[i]))
    print(subgraph1)    
subgraph_string(par,obs,W)

'''I would want capital W to be passed in the function to use ord so from 
where ever the deleted nodes are coming they should be in string...
Also when iterating over all possible W, i guess it wil be easier to have W 
originally is string'''

'''Spider doesn't recognize bitfuncs module in another file from line 9 
so using line 8'''


import pyximport
pyximport.install()
import bitfuncs
from bitfuncs import val2nice
par=[0b0,0b1,0b11,0b0,0b0,0b11000, 0b100000,0b1010100]#array of parents
obs=0b11100110 #observed nodes
W=0b1110010 #deleted nodes......this way it works even if you dont enter the nodes in alphabetical order
def subgraph(par,obs,W):
    i,new_par=0,par[:]
    while i<len(par):
        new_par[i]=par[i]&~W
        if (W & (1 << i)):
            new_par[i]=0b0
        i=i+1 
    #print(new_par)     
    return new_par   
    #print(new_par)          
subgraph(par,obs,W)   

strat=enumdags.dag_indeps(par)
print(strat)

par=[0b0,0b0,0b11,0b0,0b1110,0b10101]
obs=0b111100
X=0b100000
Y=0b1000
Z=0b100
W=0b10000

def eseparation(X,Y,Z,W,par,obs):
    d=closure.Triplet(X,Y,Z)
    new_par=subgraph(par,obs,W)
    strat=enumdags.dag_indeps(new_par)
    trips=closure.closure(strat)
    if(enumdags.implies(trips,d)==True):
        #print("True")
        return True
    else:
        #print("False")
        return False
        
eseparation(X,Y,Z,W,par,obs) 


par=[0b0,0b0,0b11,0b0,0b1110,0b10101]
obs=0b111100
X=0b100000
Y=0b1000
Z=0b10011
W=0b101000

def nondescendant_string(Z,W,obs,par):
    i=7
    j=7
    n=[]
    strat=enumdags.dag_indeps(par)
    trips=closure.closure(strat)
    summ=0
    while i>=0:
        m=0
        if(Z & (1<<i)):
            a=bitfuncs.nice2val(alph[i])
            print("*")
            while j>=0:
                if(W & (1<<j)):
                    b=bitfuncs.nice2val(alph[j])
                    c=par[j]
                    print("a")
                    d=closure.Triplet(a,b,c)
                    if(alph[i] in bitfuncs.val2nice(par[j])):
                        m=m+0
                        print("b")
                        break
                    if(alph[i] not in bitfuncs.val2nice(par[j])):
                        if(enumdags.implies(trips,d)==True):
                            m=m+0
                            print("c")
                        if(enumdags.implies(trips,d)==False):
                            m=m+1
                            print("y")
                            break
                n.append(m)
                j=j-1
               
        i=i-1        
        
    print(n)
    print(m)
    for x in n:
        summ=summ+x        
    print(summ)       
    if(summ==0):
       print("No node in Z is descendant from W")
    else:
       print("Some node in Z is descendant from W")   
    
nondescendant_string(Z,W,obs,par)   