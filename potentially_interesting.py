#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 20:55:20 2021
""""""
@author: Shashaank
"""
import itertools, multiprocessing
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags
import eseparation_functions
import bitfuncs
import math

potentially_interesting=enumdags.enumdags(5,1)    
count=0
for par in potentially_interesting:
    obs=par[1]
    parent=par[0]
    for x in enumdags.subsetsof(obs):
        for y in enumdags.subsetsof(obs & ~x):
             for z in enumdags.subsetsof(obs & ~x & ~y):
                  w= obs & ~x & ~y & ~z
                  if(len(eseparation_functions.arr(x))==1 and len(eseparation_functions.arr(y))==1):
                       if(eseparation_functions.new_theorem(x,y,z,w,parent,obs)==True):
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
      
    
                                  
                             