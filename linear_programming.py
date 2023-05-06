#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:40:38 2022

@author: sk1808
"""
import sys
import mosek
inf=0.0
def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

env = mosek.Env()
    
def linprog_mosek(objective, constraints, bounds):
    with env.Task(0, 0) as task:
        #task.set_Stream(mosek.streamtype.log, streamprinter)        
        task.appendcons(len(constraints))
        task.appendvars(len(objective))        
        for i in range(len(constraints)):
               task.putconbound(i, mosek.boundkey.up, -inf, bounds[i])
        for j in range(len(objective)):
               task.putcj(j,objective[j])                                   #constraints and variables sahre the same bound
               task.putvarbound(j, mosek.boundkey.fr, -inf , +inf) 
        for i in range(len(constraints)):
               for j in range(len(objective)):
                       if constraints[i][j]!=0:
                                task.putaij(i,j,constraints[i][j])
        task.putobjsense(mosek.objsense.minimize)
        task.optimize()
        #task.solutionsummary(mosek.streamtype.msg)
        solsta = task.getsolsta(mosek.soltype.bas)
        optimal_value=0
        if (solsta == mosek.solsta.optimal):
                xx = [0.] * len(objective)
                task.getxx(mosek.soltype.bas, xx)
                #print("Optimal solution: ")
                for i in range(len(objective)):
                        #print("x[" + str(i) + "]=" + str(xx[i]))   
                        optimal_value=optimal_value + xx[i]*objective[i]  
                return optimal_value          
        elif (solsta == mosek.solsta.dual_infeas_cer or
                solsta == mosek.solsta.prim_infeas_cer):
                #print("Primal or dual infeasibility certificate found.\n")
                return -1
        elif solsta == mosek.solsta.unknown:
                #print("Unknown solution status")
                return -1
        else:
                #print("Other solution status")
                return -1
                   
"Below is just a test"       
a=linprog_mosek([3,1,5,1],[[3,1,2,0],[2,1,3,1],[0,2,0,3]],[0,0,0]) 
print(a)           
         
            
                
