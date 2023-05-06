#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 20:43:03 2021

@author: Shashaank
"""
#Input Gdag of 4 nodes
import numpy as np
gdag=input("Enter the nodes separated by a space like: 'AB' if there is a direct edge from A to B ")
gdag_nodes= gdag.split()
print(gdag_nodes)

test=input(("Input the variables between which you want to check for condition independence; the 3rd variable ebing the conditioning one"))
test_nodes=test.split()
print(test_nodes)

#Building Ancestor Graph 
i=0
parent_graph=[]
while i<len(gdag_nodes):
    string=gdag_nodes[i]
    if (string[1]==test_nodes[0] or string[1]==test_nodes[1] or string[1]==test_nodes[2]):
        parent_graph.append(gdag_nodes[i])
    if(string not in parent_graph):
        if(string[0]==test_nodes[0] or string[0]==test_nodes[1] or string[0]==test_nodes[2]):
           parent_graph.append(gdag_nodes[i])   
    i=i+1
print(parent_graph)        

#Moralzing Ancestor Graph
moralized_graph=[]
j=0
#ask why infinite loop...when i put k=1 inside the first loop and when when k=j+1 then a repeated edge again....it should have been that when i put k=1 then the repeated edge should be there but when i put k=j+1 the repeated edge should have gone
while(j<len(parent_graph)):
     k=j+1
     while(k<len(parent_graph)):
         string1=parent_graph[j]
         string2=parent_graph[k]
         if(string1[1]==string2[1]):
            string3= string1[0]+string2[0] 
            parent_graph.append(string3)
         
         k=k+1
     j=j+1
print(parent_graph)    