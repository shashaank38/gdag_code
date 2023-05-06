#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 03:02:40 2021

@author: shashaank
"""
import math
from bitfuncs import *
from shannon_cone import *
"To remove constraints in Shanon Inequalities"   
def remove_constraints(constraints,N,allnodes,obs, shanon_cone):
    print("printing this")
    entropy_cone,i, j, k, = [], 0, 0, 0
    #for conditions in cond_indep:
        #constraints.append(conditional_mutual_information(N,conditions[0],conditions[1],conditions[2])) 
    shanon_cone_copy, constraints_copy = shanon_cone[:], constraints[:]
    #print(shanon_cone,len(shanon_cone),"the  cone before removing constraints")
    while i<len(constraints):
         if len(constraints_copy)!=0: 
             equalities ,b= [],0
             j=1
             shanon_cone=shanon_cone_copy[:]
             #print("1st if")
             if any(v!=0.0 for v in  constraints_copy[0]):
                    #print("success: inside if")
                    k=0   
                    while k<len(shanon_cone): 
                           if constraints_copy[0]!=shanon_cone[k]: 
                                 if any((allnodes & ~obs) in nonemptysubsetsof(z) and y!=0 for z, y in enumerate(constraints_copy[0])):
                                        #print("N")
                                        for index, m in enumerate(constraints_copy[0]):
                                               if m!=0 and index!=0 and (allnodes & ~obs) in nonemptysubsetsof(index):
                                                       #print(val2nice(index),"unobserved corresponding to inequalities")
                                                       entropy_cone.append(matrix_manipulate(shanon_cone[k], constraints_copy[0], index))
                                                       break;
                                 else: 
                                       # print("Y")
                                        for index, m in enumerate(constraints_copy[0]):
                                               if m!=0 and index!=0:
                                                      # print(val2nice(index),"observed corresponding to inequalities")
                                                       entropy_cone.append(matrix_manipulate(shanon_cone[k], constraints_copy[0], index))
                                                       break;
                           k=k+1 
                    shanon_cone_copy=entropy_cone[:]
                    #print(entropy_cone, len(entropy_cone))
                    entropy_cone=[] 
                    #print(entb=[[0, 0, 0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, -1, 0, 1], [0, 0, 0, -1, 0, 0, 0, 1], [-1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, -2.0]]ropy_cone, len(entropy_cone))
                    #print("no")  
                    #print(constraints_copy)
                    while j<len(constraints_copy):
                           if any((allnodes & ~obs) in nonemptysubsetsof(z) and y!=0 for z, y  in enumerate(constraints_copy[0])):
                                    for index, m in enumerate(constraints_copy[0]):
                                             if m!=0 and index!=0 and (allnodes & ~obs) in nonemptysubsetsof(index):
                                                       #print(val2nice(index), "unobserved corresponding to equalities")
                                                       equalities.append(matrix_manipulate(constraints_copy[j],constraints_copy[0], index))
                                                       break;
                           else:
                                    for index, m in enumerate(constraints_copy[0]):
                                               if m!=0 and index!=0:
                                                       #print(val2nice(index), "observed corresponding to equalities")
                                                       equalities.append(matrix_manipulate(constraints_copy[j],constraints_copy[0], index))
                                                       break;
                               
                           j=j+1  
                    constraints_copy=equalities[:]
                    #print(constraints_copy)
                    #print("yes")
                    i=i+1
             else:
                    constraints_copy.remove(constraints_copy[0])
                    i=i+1
                    continue;
         else:
              break;
         #i=i+1 
         #b=b+1
    #print(b)     
    #return (shanon_cone_copy)
    shanon_cone_copy_copy=shanon_cone_copy[:]
    y=0
    print( shanon_cone_copy,len(shanon_cone_copy_copy), "after substituting constraints in inequalities and equalities")
    "forget about this net loop. it is useless"
    
   #while y<len(shanon_cone_copy_copy):
        #print(y)
    #    z=y+1
     #   while z<len(shanon_cone_copy_copy)-1:
      #    if(y<len(shanon_cone_copy)):
       #     if  shanon_cone_copy_copy[y]==shanon_cone_copy_copy[z] or  all(v==0 for v in shanon_cone_copy_copy[y]):
        #         #print(len(shanon_cone_copy),y,z)
         #        #print("true")
          #      #print(shanon_cone_copy, "this is cone")
           #      #print(shanon_cone_copy[y])
            #     shanon_cone_copy.remove(shanon_cone_copy[y])
             #   #print(shanon_cone_copy,y)
              #   #print("true")
           # z=z+1  
          #else:
           #  break;
        #y=y+1
    print(shanon_cone_copy,len(shanon_cone_copy), "after deleting duplicate ones")
    z,y,l=0,0,0
    while z<len(shanon_cone_copy_copy):
         list_constraints=shanon_cone_copy[:]
         c=shanon_cone_copy_copy[z][:]
         b=[0 for i in range(len(list_constraints)-1)]
         list_constraints.remove(c) 
         motzkin0=list_constraints[:]
         motzkin0=[[-1* i for i in j] for j in motzkin0]
         res = linprog(c, motzkin0 , b)
         if res["fun"] * (+1)>-0.5 and res["success"]: 
             shanon_cone_copy.remove(c)
         z=z+1      
    return shanon_cone_copy;     

def eliminate(cond_indep,N,allnodes,obs):
    constraints,shanon_cone=[],Shanon_Cone(N,allnodes)
    for conditions in cond_indep:
        constraints.append(conditional_mutual_information(N,conditions[0],conditions[1],conditions[2])) 
    i=1
    while i!=0:
         new_constraints=[]
         shanon_cone_copy=remove_constraints(constraints,N,allnodes,obs,shanon_cone)
         shanon_cone_copy_negative=[[-1* y  for y in z] for z in shanon_cone_copy]
         shanon_cone=[]
         j=0
         while j<len(shanon_cone_copy):
             k=0
             while k<len(shanon_cone_copy_negative):
                 if shanon_cone_copy[j]==shanon_cone_copy_negative[k]:
                       new_constraints.append(shanon_cone_copy[j])
                 k=k+1     
             j=j+1
         if len(new_constraints)!=0:   
             constraints=new_constraints[:]
             shanon_cone=shanon_cone_copy[:]
             i=i+1
         else:
             break;
             #i=0
    return shanon_cone_copy;         
         
         
         
        

a=eliminate([[4,1,2],[8,3,4]],4,15,15)
print("in here")
print(a)
print(len(a))
def inequality(inequality,N):
    a,b=[],2**N
    for i in inequality:
        c=[]
        j=0
        while j<b:
            c.append(str(i[j])+val2nice(j))
            j=j+1
        a.append(c)    
        
    print(a)
    return a;
inequality(a,3)   
b=[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, -1.0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 1, -1, 0, 0], [-1, 1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 1, 0, -1], [0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0], [0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0]]
c=[[-1*i for i in j]for j in a]
for i in a:
    if  i not in b:
        print(i)