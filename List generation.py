#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 01:16:52 2021

@author: Shashaank
"""

import itertools, multiprocessing
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags
import eseparation_functions
import bitfuncs
import math


par=[['', '', '', 'AC', 'AB'], 'BCDE',
['', '', '', 'BC', 'AC', 'AB'], 'DEF',
['', '', '', 'BC', 'AC', 'AB'], 'ADEF',
['', '', '', 'BC', 'ACD', 'AB'], 'CDEF',
['', '', '', 'BC', 'AD', 'AB'], 'CDEF',
['', '', '', 'BC', 'AC', 'ABE', 'ABDF'], 'DEFG',
['', '', '', 'BC', 'AC', 'ABE', 'AD'], 'DEFG',
['', '', '', 'BC', 'AC', 'ABE', 'ADF'], 'DEFG',
['', '', '', 'C', 'BCD', 'ADE', 'AB'], 'CEFG',
['', '', '', 'C', 'BD', 'ACDE', 'AB'], 'CEFG',
['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'DEFG',
['', '', '', 'BC', 'AC', 'ABE', 'ABD'], 'DEFG',
['', '', '', 'BC', 'AC', 'BE', 'ADF'], 'DEFG',
['', '', '', 'BC', 'AD', 'ACE', 'ABF'], 'DEFG',
['', '', '', 'AC', 'BC', 'AE', 'ABDF'], 'DEFG',
['', '', '', 'AC', 'ABD', 'BC', 'A'], 'DEFG',
['', '', '', 'C', 'BCD', 'ACE', 'AB'], 'DEFG',
['', '', '', 'A', 'CD', 'BDE', 'BC'], 'AEFG',
['', '', '', '', 'BCD', 'ACDE', 'AB'], 'CEFG',
['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'CEFG',
['', '', '', 'BC', 'AD', 'CE', 'ABEF'], 'DEFG',
['', '', '', 'BC', 'AD', 'CE', 'ABF'], 'DEFG',
['', '', '', 'BC', 'C', 'ACDE', 'AB'], 'DEFG',
['', '', '', 'BC', 'AD', 'ACE', 'ABEF'], 'DEFG',
['', '', '', 'AC', 'BD', 'BCE', 'AEF'], 'DEFG']

print(bin(bitfuncs.nice2val("BC")))
parent=[]
obs=[]
for y in range(len(par)):
    _=[]
    __=[]
    if (y%2==0):
       for x in par[y]:
           _.append(bitfuncs.nice2val(x))
       __.append(bitfuncs.nice2val(par[y+1]))
       parent.append(_)
       obs.append(__)      
#print(parent, len(parent))
#print(obs,len(obs)) 
parent_nodes_array=[]          
for x in range(len(obs)):
    _=[]
    for y in obs[x]:
        _.append(eseparation_functions.arr(y))
        _.reverse()
        
    parent_nodes_array.append(_)
    
#print(parent_nodes_array)    
second=[]

        
        