#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:31:54 2021

@author: Shashaank
"""

parents=['','A','AB','','','DE','F','CEG']
modified=[]
i=0
#index of node to be deleted
j=5
k='F'
#Detect all the places where k='...' comes. Delete all arrows from k='...' 
#to its children
while i<len(parents):
   string1=parents[i]
   string2=''
   l=0
   while l<len(string1):
       if(string1[l]!=k):
            string2=string2+ string1[l]
           
       if(string1[l]=='C'):
           string2=string2+''
       l=l+1    
   modified.append(string2)
   i=i+1
print(modified)

#Delete the arrows from all parents of k='....' to k='....' to make it parentless
modified[j]=''
print(modified)  

a=[1,2,3]
b=[bin(a[0]), bin(a[1]), bin(a[2])]
l=0
for l in range(len(a)):
   print(a[l])
