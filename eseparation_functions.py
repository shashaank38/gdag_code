
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 18:24:31 2021

@author: Shashaank
"""
"""Program to check the conditions of the theorem in the paper for given
 nodes X, Y, Z, W, parents and observed nodes """





import pyximport # to run cython
pyximport.install()


import math
from enumdags import *
from closure import *
from bitfuncs import *





"To count the number of bits in an integer"
def countBits(number):
    if (number!=0):
       return int((math.log(number) / math.log(2)) + 1);
    else:
       return 1
 
    
 
    
 
"To get the array of bits in an integer"
def arr(A):
    count=countBits(A)    
    a=[]
    while count>=0:
        if(A & (1<<count)):
           a.append(2**count)
        count=count-1  
    return a  


#   j=j+1   


"To get the array of parents of nodes"
def parents(A,par):
    count=countBits(A)    
    parent=[]
    while count>=0:
        if(A & (1<<count)):
           parent.append(par[count])
        count=count-1  
    return parent  





" To construct the gdag with deleted nodes"
def subgraph(par,obs,W):
    i,new_par=0,par[:]
    while i<len(par): # for i in range(len(par)):
        if (W & (1 << i)):
            new_par[i]=0b0
        else:
            new_par[i]=par[i]&~W     
        i=i+1     
    return new_par   




       
"To check conditions of eseparation"    
def eseparation(X,Y,Z,W,par,obs):
    d=Triplet(X,Y,Z)
    new_par=subgraph(par,obs,W)
    strat=dag_indeps(new_par)
    trips=closure(strat)
    return implies(trips,d) 
    
    


    
"To check whether any node in Z is a descendant from W"
def nondescendant(Z,W,obs,par):
    z,w,parent,summ=arr(Z),arr(W),parents(W,par),0
    strat=dag_indeps(par)
    trips=closure(strat)
    if(Z!=0 and W!=0):
        for i in range(len(z)):
            for j in range(len(w)):
                 if(z[i] not in subsetsof(parent[j])):
                     if(implies(trips, Triplet(z[i],w[j],parent[j]))==False):
                        return False
        return True    
    else:
        return True
   
    



"To check the conditions of the NEW CORRECT (OUR'S) eseparation theorem"
def new_theorem(X,Y,Z,W,par,obs):
    strat=dag_indeps(par)
    trips=closure(strat)
    count=0
    if(X!=0 and Y!=0):
         if (eseparation(X,Y,Z,W,par,obs)==True and nondescendant(Z,W,obs,par)==True):
             for S in subsetsof(obs):
                  if(implies(trips, Triplet(X,Y,S))==True):
                      return False
             return True 
           
       
        


"To check the conditions of the OLD INCORRECT JACQUES' eseparation theorem"            
def old_theorem(X,Y,Z,W,par,obs):
    strat=dag_indeps(par)
    trips=closure(strat)
    count=0
    if(X!=0 and Y!=0):
         if (eseparation(X,Y,Z,W,par,obs)==True and nondescendant(Z,W,obs,par)==True):
             for S in subsetsof(W):
                  if(implies(trips, Triplet(X,Y,Z|S))==True):
                      return False
             return True     
 
            
 
 
   
"Finding DAGS checked and unchecked by the correct eseparation theorem" 
def potentially_interesting(n,process):
    potentially_interesting=enumdags(n,process)    
    count=0     
    for par in potentially_interesting:
       obs=par[1]
       parent=par[0]
       for x in subsetsof(obs):
            for y in subsetsof(obs & ~x):
                for z in subsetsof(obs & ~x & ~y):
                    w= obs & ~x & ~y & ~z
                    if(len(arr(x))==1 and len(arr(y))==1):
                         if(new_theorem(x,y,z,w,parent,obs)==True):
                             print([[val2nice(i) for i in parent], val2nice(obs)],"   is interesting")
                             print("The required (X,Y,Z,W) are:","(", val2nice(x), ",",val2nice(y),",",val2nice(z),",",val2nice(w),")")
                             count=count+1
                             break
                else:
                   continue
                break
            else:
               continue
            break
    else:           
        print([[val2nice(i) for i in parent], val2nice(obs)],"    could not be checked")
    print("Number of interesting GDAGS that could be checked is: ",count)  





"Calling for final output"
#n=int(input("Enter the number of nodes, the dags of which you want to check esparation method on; parallel processes=1"))   
#potentially_interesting(n, 1)        
         
