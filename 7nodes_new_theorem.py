#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 13:20:53 2021

@author: Shashaank
"""

import itertools, multiprocessing
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags
import eseparation_functions
import bitfuncs
import math


par= [['', '', 'AB', 'AC'], 'BCD',
['', '', 'AB', 'BC', 'AC'], 'CDE',
['', '', '', 'AC', 'AB'], 'BCDE',
['', '', 'AB', '', 'BD', 'ACE'], 'CDEF',
['', '', '', 'BC', 'AD', 'ABE'], 'CDEF',
['', '', 'A', 'BC', 'CD', 'BE'], 'ADEF',
['', '', 'AB', '', 'BCD', 'AE'], 'CDEF',
['', '', '', 'BC', 'AC', 'AB'], 'DEF',
['', '', 'AB', 'AC', 'BD', 'AD'], 'CDEF',
['', '', 'AB', 'AC', 'BD', 'A'], 'CDEF',
['', '', 'B', 'ABC', 'BD', 'AE'], 'CDEF',
['', '', 'B', 'ABC', 'BCD', 'AE'], 'BDEF',
['', '', 'AB', 'B', 'BCD', 'AE'], 'CDEF',
['', '', '', 'BC', 'AC', 'AB'], 'ADEF',
['', '', 'B', 'ABC', 'CD', 'AE'], 'BDEF',
['', '', '', 'BC', 'ACD', 'AB'], 'CDEF',
['', '', 'AB', '', 'BCD', 'ACE'], 'CDEF',
['', '', '', 'ABC', 'BCD', 'AE'], 'BDEF',
['', '', 'B', 'AC', 'BCD', 'AE'], 'BDEF',
['', '', '', 'BC', 'AD', 'AB'], 'CDEF',
['', '', 'B', 'ABC', 'BCD', 'AE'], 'CDEF',
['', '', 'AB', 'C', 'CD', 'BCDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'CD', 'BDE', 'ACF'], 'CEFG',
['', '', '', 'BC', 'AC', 'ABE', 'ABDF'], 'DEFG',
['', '', 'B', 'BC', 'ABCD', 'BDE', 'AF'], 'CEFG',
['', '', '', 'BC', 'AC', 'BDE', 'AF'], 'DEFG',
['', '', '', 'BC', 'ACD', 'BE', 'AD'], 'DEFG',
['', '', '', 'BC', 'AC', 'ABE', 'AD'], 'DEFG',
['', '', 'B', 'BC', 'ABCD', 'CDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'BCD', 'ACE', 'D'], 'CEFG',
['', '', '', 'BC', 'AC', 'ABE', 'ADF'], 'DEFG',
['', '', '', 'C', 'BCD', 'ADE', 'AB'], 'CEFG',
['', '', '', 'BC', 'AD', 'CE', 'AB'], 'DEFG',
['', '', '', 'C', 'BD', 'ACDE', 'AB'], 'CEFG',
['', '', 'B', 'BC', 'ACD', 'BCDE', 'AF'], 'CEFG',
['', '', '', 'BC', 'ACD', 'BE', 'AF'], 'DEFG',
['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'DEFG',
['', '', 'AB', 'BC', '', 'ABDE', 'ACDF'], 'EFG',
['', '', 'AB', 'C', 'BCD', 'DE', 'AE'], 'CEFG',
['', '', 'AB', 'C', 'D', 'BDE', 'ACF'], 'CEFG',
['', '', 'AB', 'C', 'BCD', 'DE', 'ACE'], 'CEFG',
['', '', '', 'BC', 'AC', 'ABE', 'ABD'], 'DEFG',
['', '', '', 'BC', 'AC', 'BE', 'ADF'], 'DEFG',
['', '', '', 'BC', 'AD', 'ACE', 'ABF'], 'DEFG',
['', '', 'B', 'BC', 'ABD', 'DE', 'AF'], 'CEFG',
['', '', 'B', 'BC', 'ABD', 'CDE', 'AF'], 'CEFG',
['', '', '', 'BC', 'ACD', 'BDE', 'AD'], 'DEFG',
['', '', '', 'AC', 'BD', 'BCE', 'AE'], 'DEFG',
['', '', '', 'ABC', 'CD', 'BE', 'AE'], 'DEFG',
['', '', '', 'AC', 'BC', 'AE', 'ABDF'], 'DEFG',
['', '', 'B', 'C', 'ABCD', 'BDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'BD', 'AE', 'D'], 'CEFG',
['', '', 'B', 'C', 'ABD', 'BDE', 'AF'], 'CEFG',
['', '', 'B', 'BC', 'AD', 'BCDE', 'AF'], 'CEFG',
['', '', '', 'AC', 'ABD', 'BC', 'A'], 'DEFG',
['', '', 'AB', 'BC', '', 'ACDE', 'ABDF'], 'EFG',
['', '', 'AB', 'C', 'BD', 'AE', 'CD'], 'CEFG',
['', '', 'AB', 'C', 'CD', 'BDE', 'AF'], 'CEFG',
['', '', 'AB', '', 'CD', 'BCDE', 'ACF'], 'CEFG',
['', '', '', 'AC', 'BD', 'BCDE', 'AF'], 'DEFG',
['', '', '', 'C', 'BCD', 'ACE', 'AB'], 'DEFG',
['', '', 'AB', 'C', 'BCD', 'ACE', 'CD'], 'CEFG',
['', '', 'B', 'C', 'ABCD', 'BCDE', 'AF'], 'CEFG',
['', '', '', 'BC', 'AC', 'BDE', 'AD'], 'DEFG',
['', '', 'B', 'BC', 'ABD', 'BDE', 'AF'], 'CEFG',
['', '', '', 'A', 'CD', 'BDE', 'BC'], 'AEFG',
['', '', 'AB', 'C', 'BCD', 'AE', 'CD'], 'CEFG',
['', '', '', '', 'BCD', 'ACDE', 'AB'], 'CEFG',
['', '', '', 'BC', 'AC', 'BE', 'AD'], 'DEFG',
['', '', 'AB', 'C', 'BD', 'ACE', 'CD'], 'CEFG',
['', '', 'B', 'BC', 'AD', 'BDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'BCD', 'AE', 'D'], 'CEFG',
['', '', '', 'BC', 'ACD', 'BDE', 'AF'], 'DEFG',
['', '', 'AB', 'C', 'D', 'BDE', 'AF'], 'CEFG',
['', '', 'AB', '', 'BCD', 'CDE', 'ACE'], 'CEFG',
['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'CEFG',
['', '', 'AB', 'C', 'D', 'BCDE', 'AF'], 'CEFG',
['', '', 'B', 'C', 'ABD', 'BCDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'BCD', 'CDE', 'ACE'], 'CEFG',
['', '', 'AB', 'C', 'D', 'BCDE', 'ACF'], 'CEFG',
['', '', '', 'BC', 'AD', 'CE', 'ABEF'], 'DEFG',
['', '', 'B', 'BC', 'ABCD', 'DE', 'AF'], 'CEFG',
['', '', '', 'BC', 'AD', 'CE', 'ABF'], 'DEFG',
['', '', 'AB', 'C', 'BD', 'ACE', 'D'], 'CEFG',
['', '', '', 'BC', 'C', 'ACDE', 'AB'], 'DEFG',
['', '', 'AB', '', 'BD', 'CDE', 'ACE'], 'CEFG',
['', '', 'AB', 'C', 'BD', 'CDE', 'AE'], 'CEFG',
['', '', 'B', 'BC', 'ACD', 'BDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'CD', 'BCDE', 'ACF'], 'CEFG',
['', '', 'AB', 'C', 'BD', 'DE', 'ACE'], 'CEFG',
['', '', '', 'AC', 'BD', 'BCE', 'AF'], 'DEFG',
['', '', '', 'BC', 'AD', 'ACE', 'ABEF'], 'DEFG',
['', '', 'B', 'BC', 'ABD', 'BCDE', 'AF'], 'CEFG',
['', '', 'AB', 'C', 'BD', 'CDE', 'ACE'], 'CEFG',
['', '', '', 'AC', 'BD', 'BCE', 'AEF'], 'DEFG',
['', '', 'AB', '', 'CD', 'BDE', 'ACF'], 'CEFG',
['', '', 'AB', 'C', 'BCD', 'CDE', 'AE'], 'CEFG',
['', '', 'B', 'BC', 'ABCD', 'BCDE', 'AF'], 'CEFG',
['', '', '', 'AC', 'BCD', 'BDE', 'AE'], 'DEFG',
['', '', 'AB', 'C', 'BD', 'DE', 'AE'], 'CEFG']


count=0
for y in range(len(par)):
    if (y%2==0):
       parent=[]
       obs=0
       for x in par[y]:
           parent.append(bitfuncs.nice2val(x))
       obs=bitfuncs.nice2val(par[y+1])
       for x in enumdags.subsetsof(obs):
        for y in enumdags.subsetsof(obs & ~x):
             for z in enumdags.subsetsof(obs & ~x & ~y):
                  w= obs & ~x & ~y & ~z
                  if(len(eseparation_functions.arr(x))==1 and len(eseparation_functions.arr(y))==1):
                       if(eseparation_functions.new_theorem(x,y,z,w,parent,obs)==True ):
                           print([[bitfuncs.val2nice(i) for i in parent], bitfuncs.val2nice(obs)],"is interesting")
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
        print([[bitfuncs.val2nice(i) for i in parent], bitfuncs.val2nice(obs)])
print("Number of interesting GDAGS that could be checked is: ",count)  
print("Number of not checked GDAGS is: ",len(par)/2-count)  
           
      