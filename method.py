#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 23:54:00 2021

@author: Shashaank
"""
#parents=['','A','A','BC','BD','CE','DEF']
check=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
string1='EG'
string2='AC'
i=0
while i<len(string1):
    string3=string1[i]
    j=0
    while j<len(check):
        if(check[j]==string3):
           print(j) 
        j=j+1                #Why if I delete this line then it prints nothing just keeps on looping but if i delete the below line of i=i+1 then it keeps on printing instead
    print(j)       
    i=i+1     
    #while j<len(string2):
a='P'        
i=ord(a)-65
print(i)
i=0
j=0
while i<4:
    if(i==3):
        print(i)
        break
    else:
        print(i)
    j=i+1 
    print(j)
    i=i+1        
    
from bitfuncs import val2nice
from bitfuncs import nice2val
from enumdags import dag_indeps
from closure import closure
from enumdags import implies
from closure import Triplet
a=0b1010
c=val2nice(a)
print(c)    

par=[0b0,0b1,0b10,0b100,0b1000]
k='E'
c='A'
def non_descendent():
    strat=dag_indeps(par)
    trips=closure(strat)
    w=k
    z=c
    i=0
    m=0
    n=[]
    o=0
    while i<len(z):
        string1=z[i]
        d=nice2val(string1)
        j=0
        while j<len(w):
            string2=w[j]
            e=nice2val(string2)
            position=ord(string2)-65
            string3=val2nice(par[position])
            f=nice2val(string3)
            g=Triplet(d,e,f)
            if(implies(trips,g)==True):
                m=0
            if(implies(trips,g)==False):    
                m=1
                break
            j=j+1
        n[i]=m    
        i=i+1
    for x in n:
        if(x==0):
            o=0
        else:
            o=1
    if(o==0):
       print("No node in Z is descedent from W")
    else:
       print("Some node in Z is descedent from W")         
            
            
non_descendent()                