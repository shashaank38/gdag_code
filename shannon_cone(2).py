#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 17:22:53 2021

@author: Shashaank
"""

import pyximport # to run cython
pyximport.install()
import time
import timeit
import math
from enumdags import *
from closure import *
from bitfuncs import *
from eseparation_functions import *
from itertools import chain, combinations
from scipy.optimize import linprog


"To generate the power set of a set"
def powerset(iterable):
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))


"To create a list of gdag variables in the notation of the earlier code"           
def dag(N):
  i=0 
  dag=[]
  while i<N:
      dag.append(2**i)
      i=i+1
  return dag 
 

"""Generator of the bitstrings with a single 1 in the place 
of each 1 in bitstring left"""
def bitsof_inbits(left):
    x=0
    while left:
        x = left & ~(left - 1)
        left &= ~x
        yield x
        

"To turn off the $ K^{th} $ bit in an integer $ N $"
def turnOffK(n,k):
    if (k <= 0):
        return n
    return (n & ~(1 << (k - 1))) 


def zero_list(list_element):
    a=0
    for i in list_element:
        if i!=0:
            a=a+1
    return a        
  


"""To generate $ N $ shannon inequalties following from positivity of 
conditional entropy: $ H(X_i|X_{N-i}) >= 0 $ given the total number of 
nodes $ N $"""

def conditional_entropy_vectors(N):
    dag_var=dag(N)
    all_vec=[[0]*(2**N) for i in range(N)] 
    k=0
    for i in dag_var:
        l,m=0,0
        for j in range(len(dag_var)):
            if dag_var[j]!=i:
                l=l|(dag_var[j]) 
            m=m|(dag_var[j])  
        all_vec[k][l]=-1
        all_vec[k][m]=1
        k=k+1      
    return all_vec        
           
    

"""To generate $ C^N_2*(2^{n-2}) $ shannon inequalties following from positivity of 
condition mutual information: $ I(X_i; X_j| X_K) >= 0 $, where $ K $ is 
subset of $ N-{i,j} $ given the total number of nodes $ N $"""

def conditional_mutual_information_vectors(N):
  if(N>1): 
    dag_var=dag(N)
    all_vec=[[0]*(2**N) for i in range(int((math.factorial(N)/(math.factorial(N-2)*2)*2**(N-2))))]
    s,x,j=0,0,0
    for i in dag_var:
       s=s|i   
    while j<len(dag_var):
         k=j+1
         while k<len(dag_var):
              for l in subsetsof(s & ~(dag_var[j]) & ~(dag_var[k])):
                    m,n,o,p=0,0,0,0 
                    m=m|dag_var[j]|dag_var[k]|l
                    n=n|l
                    o=o|dag_var[k]|l
                    p=p|dag_var[j]|l
                    all_vec[x][m]=-1
                    all_vec[x][n]=-1
                    all_vec[x][o]=1
                    all_vec[x][p]=1
                    x=x+1
              k=k+1
         j=j+1        
    return all_vec
  else:
      all_vec=[]
      return all_vec


"""Note that in the above 2 functions the 1st component of each vector
 denotes the coefficient of $ H(\emptyset) $ in the expansion. Since 
 $ H(\emptyset) $ hits in particular expansions, the corresponding 
 entropy vectors have a $\pm 1$ in the first component. It doesn't mean 
 that H(Null) $ \neq $ 0"""
 
 
"""The above 2 functions return all the Shanon type inequalities for a given 
$ N $ Of course you need to append both the lists together""" 



"""To get A conditional entropy vector's components 
$ ( H(X|Y) = H(XY) - H(Y) ) $ """
def conditional_entropy(N,X,Y):
   vec = [0] * (1<<N)
   vec[X|Y]=1
   vec[Y]=-1
   return vec    



"""To get all $ N $ Shanon type entropy vectors following from positivity of 
conditional entropy"""
def Conditional_Entropies(N, allnodes):
    #allnodes = (1<<N)-1;
    return [conditional_entropy(N,i,allnodes & ~i) for i in bitsof_inbits(allnodes)];




"""To get A conditional mutual information entropy vector's components
$ ( I(X;Y|Z) = H(X,Z) + H(Y,Z) - H(Z) - H(X,Y,Z) ) $"""
def conditional_mutual_information(N,X,Y,Z):
    vec= [0] * (1<<N)
    vec[X|Z]=1
    vec[Y|Z]=1
    vec[Z]=-1
    vec[X|Y|Z]=-1
    return vec;



"""To get all the $ C^N_2*(2^{n-2}) $ Shanon type vectors following from 
positivity of conditional mutual information"""
def Conditional_Mutual_Informations(N, allnodes):
    #allnodes=(1<<N)-1; 
    allnodes1=allnodes
    vec=[]
    for i in bitsof_inbits(allnodes):
        for j in bitsof_inbits(allnodes & ~i):
            for k in subsetsof(allnodes1 & ~i & ~j):
                vec.append(conditional_mutual_information(N,i,j,k));
        allnodes=allnodes & ~i; #working here but not below
    return vec;   
     


"To get the full Shanon cone for number of nodes N"   
def Shanon_Cone(N, allnodes):
    list1=Conditional_Entropies(N, allnodes)
    list2=Conditional_Mutual_Informations(N, allnodes)
    shanon_cone = [y for x in [list1, list2] for y in x] #cooler way of doing list1 + list2
    return shanon_cone
#print(Shanon_Cone(3,7))
#print("hello")
    

"To remove variable at postion = position which is in matrix2 from matrix1"                
def matrix_manipulate(matrix1, matrix2, position): #variable=position
    l = len(matrix1) - len(matrix2)
    if l>0:
       for i in range(len(matrix2), len(matrix2) + abs(l)):
          matrix2.append(0)      
    if l<0:      
       for i in range(len(matrix1), len(matrix1) + abs(l)):
          matrix1.append(0)       
    k=matrix1[position]/matrix2[position]
    if k>0.0: 
       matrix=[matrix1[j] - (abs(k) * matrix2[j]) for j in range(len(matrix1))]
       return matrix;
    if k<0.0:
       matrix=[matrix1[j] + (abs(k) * matrix2[j]) for j in range(len(matrix1))]                      
       return matrix;
    if k==0.0:
       return matrix1; #Working Perfectly
   
    


"To remove constraints in Shanon Inequalities"   
def remove_constraints(cond_indep,N,allnodes):
    entropy_cone, constraints, shanon_cone, i, j, k, = [], [], Shanon_Cone(N,allnodes), 0, 0, 0
    for conditions in cond_indep:
        constraints.append(conditional_mutual_information(N,conditions[0],conditions[1],conditions[2])) 
    shanon_cone_copy, constraints_copy = shanon_cone[:], constraints[:]
    print(shanon_cone,len(shanon_cone),"the I cone before rmoving constraints")
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
                                 for index, m in enumerate(constraints_copy[0]):
                                         if m!=0 and index!=0:
                                                 entropy_cone.append(matrix_manipulate(shanon_cone[k], constraints_copy[0], index))
                                                 break;
                           k=k+1 
                    shanon_cone_copy=entropy_cone[:]
                    #print(entropy_cone, len(entropy_cone))
                    entropy_cone=[] 
                    #print(entropy_cone, len(entropy_cone))
                    #print("no")  
                    #print(constraints_copy)
                    while j<len(constraints_copy):
                             for index, n in enumerate(constraints_copy[0]):
                                     if n!=0 and index!=0:
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
    while y<len(shanon_cone_copy_copy):
        #print(y)
         z=y+1
         while z<len(shanon_cone_copy_copy)-1:
           if(y<len(shanon_cone_copy)):
             if  shanon_cone_copy_copy[y]==shanon_cone_copy_copy[z] or  all(v==0 for v in shanon_cone_copy_copy[y]):
                  print(len(shanon_cone_copy),y,z)
                  #print("true")
                 #print(shanon_cone_copy, "this is cone")
                  print(shanon_cone_copy[y])
                  shanon_cone_copy.remove(shanon_cone_copy[y])
                 #print(shanon_cone_copy,y)
                  #print("true")
             z=z+1  
           else:
              break;
         y=y+1
   #print(shanon_cone_copy)
    return shanon_cone_copy;     

#a=remove_constraints([[1,6,0],[1,2,0],[1,2,4],[1,4,0],[1,4,2]],3,7)
#print("in here")
#print(a)
#print(len(a))




"To get the list of all observed conditional independences in the GDAG"
def Iobservable_indeps(par,obs):
    res = []
    add=[]
    for i in nonemptysubsetsof(obs):
        for j in nonemptysubsetsof(obs & ~i):
            for k in subsetsof(obs & ~i & ~j) :
                ijk = Triplet(i,j,k)
                if implies(par, ijk):
                    res.append(ijk)
                    add.append([i,j,k])
    return add #it was "return tuple(res)" here



"To get the list of all Markov conditions in the GDAG"
def Markov_cond(par,obs,N):
     new_par,add=[],[]
     k=0
     strat=dag_indeps(par)
     #print(strat)
     trips=closure(strat)
     allnodes=(1<<N)-1; 
     for i in bitsof_inbits(allnodes):
         new_par.append(parents(i,par))
         for j in nonemptysubsetsof(allnodes & ~i & ~new_par[k][0]):
             #if(i not in subsetsof(new_par[k][0])): #is this line or the above third term after 2nd & really needed at all 
                 if(implies(trips, Triplet(i,j,new_par[k][0]))==True):
                     add.append([i,j,new_par[k][0]])
                     break;
         k=k+1            #Z is independent of previous nodes given parents of Z
     return add 
 
 

"To generate the inequalities corresponding to $ \mathcal{I} $"
def I_cone(par,obs,N):    #to correct
    strat=dag_indeps(par)
    trips=closure(strat)
    cond_indep=Iobservable_indeps(trips, obs)
    print(cond_indep, "conditional independences for I cone")
    return remove_constraints(cond_indep,N,obs)       
    
    
 
    
"To generate the inequalities corresponding to $ \mathcal{C} $"    
def C_cone(par,allnodes,obs,N):
    cond_indep=Markov_cond(par,obs,N)
    print(cond_indep, "for C cone")
    to_remove_unobserved=remove_constraints(cond_indep,N,allnodes)
   #print(to_remove_unobserved,len(to_remove_unobserved))
   #print("This is the C_cone after removing constraints....Looks good because its has some elements in decimals now.....and the remove constraint was atleast doing the calculations correctly" )
    #j=0
    #print(val2nice(allnodes & ~obs))
    to_remove_unobserved_copy=to_remove_unobserved[:]
    z,y,l=0,0,0
    while z<len(to_remove_unobserved):
        list_constraints=to_remove_unobserved_copy[:]
        c=to_remove_unobserved[z][:]
        b=[0 for i in range(len(to_remove_unobserved_copy)-1)]
        list_constraints.remove(c) 
        motzkin0=list_constraints[:]
        motzkin0=[[-1* i for i in j] for j in motzkin0]
        print("linear program number=",l)
        start=time.time()
        res = linprog(c, motzkin0 , b)
        print("time taken=", time.time()-start)
        l=l+1
       #print(res['fun'], res['success'], res['status'], res['message'])
        if res["fun"] * (+1)>-0.5 and res["success"]: 
            print(c)
            to_remove_unobserved_copy.remove(c)
       #print(z)    
        z=z+1    
    to_remove_unobserved=to_remove_unobserved_copy[:]    
   #print(to_remove_unobserved, len(to_remove_unobserved))
   #print("This is the C_cone after removing redundant inequalities")     
    
    for i in bitsof_inbits(allnodes & ~obs):
        j=i
        while j<= (2**N):
            if i in nonemptysubsetsof(j):
                unobserved_removed=[]
                unobserved_removed.append(Fourier_Motzkin(to_remove_unobserved, j))
                to_remove_unobserved=unobserved_removed[0][:]
            j=j+1    
     
    #print(to_remove_unobserved, len(to_remove_unobserved))
    #print("no")
    return to_remove_unobserved;  
            




"""Fourier-Motzkin elimination eliminating variable at matrix[position] 
from the list of all inequalities in list matrix"""    
def Fourier_Motzkin(matrix, position):
    i, motzkin, j = 0, [], 0
    while i<len(matrix):
        j=i+1
        while j<len(matrix):
            if matrix[i][position] * matrix[j][position]<0:
                motzkin.append(matrix_manipulate(matrix[i], matrix[j], position))
            j=j+1
        if matrix[i][position]==0:
            motzkin.append(matrix[i])    
        i=i+1 
    motzkin1, k,l = motzkin[:], 0,1
    while k<len(motzkin):
        start=0.0
        c=[i * 1 for i in motzkin1[k]]
        #d=[i * -1 for i in motzkin1[k]]
        motzkin1.pop(k)
        if len(motzkin1)!=0:
             b=[0 for i in range(len(motzkin1))]
             #print(c,motzkin1,b)
             #b.append(1)
             motzkin1=[[-1* i for i in j] for j in motzkin1]
             #motzkin1.append(d)
             #print(len(motzkin1), "is the size of the current constraint matri to linprog")
             print("linear program number=",l)
             start=time.time()
             res = linprog(c, motzkin1, b, method='revised simplex')
             print("time taken=", time.time()-start)
             l=l+1
             print(motzkin1)
            #print(res['fun'], res['success'], res['status'], res['message'])
             if res["fun"] * (+1)>-0.5 and res["success"]:
                 print(motzkin[k])
                 motzkin.pop(k)
                 motzkin1=motzkin[:]
                 
                #print(res['fun'], res['success'])
             else:
                 k=k+1
                 motzkin1=motzkin[:]
             for z in motzkin1:
                 if all(v==0.0 for v in z):
                     motzkin.remove(z)
                     motzkin1=motzkin[:]
        else:
             k=k+1;
    return motzkin;  

    
   
    
#Working - tests Fourier_Motzkin and matrix_manipulate functions
#a=[[0,0,1,5],[0,1,0,6],[1,0,0,4],[0,1,1,2]]
#c=0 #shuffle the above values according to the comment in green below to see the effect
#d=Fourier_Motzkin(a,c)
#print(d)
#print(len(d))
    
    

"""When the number of inequalities are quite large a VERY BIG ADVANTAGE CAN be
got from eliminating that variable in each step for which the  total number of pairs
of inequalities in which the variable occurs with opposite signs is the least.
another advantage is in eliminating that variable for which the maximum number
of inequalities has its corresponding coefficient 0"""    


                     
#Testing Markov_cond- Working (with strat==add for all the dgas I tested)
#print(Markov_cond([0b0,0b1,0b10,0b110,0b1101,0b11101,0b110111],0b1011011,7))


"To list inequalities in $\mathcal{C}$ but not in $mathcal\{I}"
def gdag_cone(par,allnodes,obs,N):
    list_C_cone=C_cone(par,allnodes,obs,N)
    list_I_cone=I_cone(par,obs,N) 
    print(list_C_cone, len(list_C_cone))
    print("This is the C_cone after Fourier Motzkin elimination of all unobserved variables' occurences...and it looks correct because we tested Fourier Motzkin individually and it gave the correct answers")
    print(list_I_cone, len(list_I_cone))
    print("This is the I_cone after removing constraints...again correct because no number is in decimals showing it never entered remove_constraints since there are no observed conditional independences")
    interesting_inequalities, non_interesting=[], []
    for i in list_C_cone:
        if i not in list_I_cone:
                interesting_inequalities.append(i)
        else:
            non_interesting.append(i)
   #print(non_interesting, len(non_interesting))
   #print("This is the list  of non interesting inequalities that are in C and I") 
    return  interesting_inequalities      
            
a=gdag_cone([0b0,0b1,0b0],0b111,0b111,3)
print(a,len(a)) 
print("This is the list of inequalities that are in C_cone but not in I_cone....They are 16...and there are no occurences of unobserved vraibles") 
c=0  
        
#for i in a:
    #if i[1]==0 and i[3]==0 and i[5]==0 and i[7]==0 and i[9]==0 and i[11]==0 and i[13]==0 and i[15]==0:
        #print("coefficients of A, AB, AC, ABC, AD, ABD, ACD, ABCD all are 0")
        #c=c+1
   
    #else:
        #print("in", i, "coefficients of A, AB, AC, ABC, AD, ABD, ACD, ABCD all are not 0")
#print(c)
