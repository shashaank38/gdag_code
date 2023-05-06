#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:34:21 2021

@author: Shashaank
"""

import itertools, multiprocessing
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags
import eseparation_functions
import bitfuncs
import math

a=enumdags.find_candidates(5, False)
par=[]
obser=[]
for i in a:
   par.append(list(i[1]))
   obser.append(i[2])  
count=0   
for i in range(len(par)):
    parent=par[i]
    obs=obser[i]
    for x in enumdags.subsetsof(obs):
        for y in enumdags.subsetsof(obs & ~x):
             for z in enumdags.subsetsof(obs & ~x & ~y):
                  w= obs & ~x & ~y & ~z
                  if(x|y|z|w==obs and x!=0 and y!=0):
                       if(eseparation_functions.old_theorem(x,y,z,w,parent,obs)==True):
                           print([[bitfuncs.val2nice(i) for i in parent], bitfuncs.val2nice(obs)],"   is interesting")
                           print("The required (X,Y,Z,W) are:","(", bitfuncs.val2nice(x), ",",bitfuncs.val2nice(y),",",bitfuncs.val2nice(z),",",bitfuncs.val2nice(w),")")
                           count=count+1
                           break
             else:
                 continue
             break
        else:
            continue
        break
    else:           
        print([[bitfuncs.val2nice(i) for i in parent], bitfuncs.val2nice(obs)],"    could not be checked")
print("Number of interesting GDAGS that could be checked is: ",count)                
     
    
                                  
                             
   

  