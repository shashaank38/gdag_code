#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 17:02:31 2022

@author: shashaank
"""

"""Classical Cutting Stock Problem solution for 12 total rolls where client demands 14 types of 
rolls of different lenghts.....The answer returns that only 9 out of 12 roles need to be used"""

from mip import Model, xsum, BINARY, INTEGER

n = 12  # maximum number of bars
L = 67200 #total bar lenghth  bar length should be greater than sum over (the length of each cut * number of rolls)
m = 14  # number of requests by client
w = [1380, 1520, 1560, 1710,1820, 1880, 1930, 2000, 2050, 2100, 2140, 	2150, 2200, 8000]  # size of each item
b = [22, 25, 12, 14, 18, 18, 20, 10, 12, 14, 16, 18, 20, 18]  # demand for each item

# creating the model
model = Model()
x = {(i, j): model.add_var(obj=0, var_type=INTEGER, name="x[%d,%d]" % (i, j))
     for i in range(m) for j in range(n)}
y = {j: model.add_var(obj=1, var_type=BINARY, name="y[%d]" % j)
     for j in range(n)}

# constraints
for i in range(m):
    model.add_constr(xsum(x[i, j] for j in range(n)) >= b[i])
for j in range(n):
    model.add_constr(xsum(w[i] * x[i, j] for i in range(m)) <= L * y[j])

# additional constraints to reduce symmetry
for j in range(1, n):
    model.add_constr(y[j - 1] >= y[j])

# optimizing the model
model.optimize()

# printing the solution
print('')
print('Objective value: {model.objective_value:.3}'.format(**locals()))
print('Solution: ', end='')
for v in model.vars:
    if v.x > 1e-5:
        print('{v.name} = {v.x}'.format(**locals()))
        print('          ', end='')