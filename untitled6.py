#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 12:37:20 2023

@author: shashaank
"""

class Addition:
    first = 0
    second = 0
    answer = 0

	# parameterized constructor
    def __init__(self, f, s):
        self.first = f
        self.second = s

    def display(self):
        print("First number = ", self.first)
        print("Second number = ", self.second)        
        print("answer=", self.answer)
        
    def calculate(self):
        self.answer = self.first + self.second


# creating object of the class
# this will invoke parameterized constructor
obj1 = Addition(1000, 2000)

# creating second object of same class
obj2 = Addition(10, 20)

# perform Addition on obj1
obj1.calculate()

# perform Addition on obj2
obj2.calculate()

# display result of obj1
obj1.display()

# display result of obj2
obj2.display()
def abc():
    a=5
    qwe(a)
def qwe(a): 
    print(a)  
    return 1