#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:45:08 2023

@author: shashaank
"""

from random import *
from linear_programming import linprog_mosek
from github_FM_advantage import *


for j in range(1):
        print(j+1,"th Polytope")
        print("Calculations for", (j+1),"th polytope begin")
        i=0
        inequalities=[]
        while i<=10:
          inequality=[randint(-5,5) for i in range(10)]
          b=[0 for i in range(len(inequalities))]
          res = linprog_mosek(inequality, inequalities, b)
          if res < -0.5:
            inequalities.append(inequality)  
          i=i+1
        print((inequalities))

