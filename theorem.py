#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 13:59:17 2021

@author: Shashaank
"""

from subgraph_bitwise import *
from non_descedant_bitwise import *
from enumdags import *
from closure import *
from bitfuncs import *


def theorem(X,Y,Z,W,par,obs):
    strat=enumdags.dag_indeps(par)
    trips=closure.closure(strat)
    if(X|Y|Z|W==obs):
       if (eseparation(X,Y,Z,W,par,obs)==True and nondescendant(Z,W,obs,par)==True):
           if(implies(trips, Triplet(X,Y,Z))==False and implies(trips, Trips(X,Y,(Z|subsetsof(W))))==False):
               #print(par)
               return True
           
              
        