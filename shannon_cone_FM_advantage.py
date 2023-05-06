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
from linear_programming import linprog_mosek
from itertools import starmap
from multiprocessing import Pool

start_time = time.time()

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

vec=conditional_entropy(6,0b011000,0b100)
print(vec)
print(val2nice(40))


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
    

def matrix_round(matrix):
  round_matrix=[]  
  for i in matrix:
        if abs(i)<0.000000001:
            round_matrix.append(0)
        else:
            round_matrix.append(i)
  return round_matrix          
    

"To remove variable at postion = position which is in matrix2 from matrix1"                
def matrix_manipulate(matrix1, matrix2, position): #variable=position
    l = len(matrix1) - len(matrix2)
    if l>0:
       for i in range(len(matrix2), len(matrix2) + abs(l)):
          matrix2.append(0)      
    if l<0:      
       for i in range(len(matrix1), len(matrix1) + abs(l)):
          matrix1.append(0) 
    if abs(matrix2[position])<0.01:
         print("Error", matrix2[position], matrix1[position],matrix1[position]/matrix2[position])  
        # print(matrix1, "yyyy")
        # print(matrix2, "nnnn")
        # k=matrix1[position]/matrix2[position]
        # print([matrix1[j] - (abs(k) * matrix2[j]) for j in range(len(matrix1))])
    k=matrix1[position]/matrix2[position]
    if k>0.0: 
        matrix=[matrix1[j] - (abs(k) * matrix2[j]) for j in range(len(matrix1))]
        return matrix_round(matrix);
    if k<0.0:
        matrix=[matrix1[j] + (abs(k) * matrix2[j]) for j in range(len(matrix1))]                      
        return matrix_round(matrix);
    if k==0.0:
        return matrix_round(matrix1); #Working Perfectly
   
    



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
                    #print(entropy_cone, len(entropy_cone))
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
    print(shanon_cone_copy,len(shanon_cone_copy),"after linprog")     
    return shanon_cone_copy; 




    
"To eliminate redundant inequalities in the  Shannon Cone"
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
         j,c=0,0
         while j<len(shanon_cone_copy):
             k=j+1
             while k<len(shanon_cone_copy_negative):
                 if shanon_cone_copy[j]==shanon_cone_copy_negative[k]:
                       new_constraints.append(shanon_cone_copy[j])
                       c=j
                       print(j)
                 k=k+1     
             j=j+1
         print(new_constraints,c)    
         if len(new_constraints)!=0:   
             constraints=new_constraints[:]
             shanon_cone=shanon_cone_copy[:]
             i=i+1
             print(i)
         else:
             i=0
             print(i)
    return shanon_cone_copy;         



        

def remove_redundancies(cond_indep,N,allnodes,obs):
    constraints,shanon_cone=[],Shanon_Cone(N,allnodes)
    negative_shanon_cone=[[-1*i for i in j]for j in shanon_cone]
    for conditions in cond_indep:
        negative_shanon_cone.append(conditional_mutual_information(N,conditions[0],conditions[1],conditions[2])) 
    total_list=negative_shanon_cone[:]
    negative_shanon_cone_copy=negative_shanon_cone[:]
    shanon_cone_copy=shanon_cone
    i=0
    for inequalities in negative_shanon_cone:
        a=inequalities
        a=[-1*i for i in a]
        b=total_list[:]
        b.pop(i)
        c=[0*i for i in range(len(b))]
        res=linprog_mosek(a,b,c)
        if res > -0.5: 
            total_list.remove(inequalities)
            #total_list.pop(i)
        else:
            i=i+1            
    total_list=[[-1* i for i in j] for j in total_list]       
    return total_list       
            
   
            
        
    

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
    return remove_redundancies(cond_indep,N,obs,obs)    



def advantange1(inequalities_list, variables_to_eliminate):
    var=[]
    count=[]
    for i in variables_to_eliminate:
        j=0
        c=0
        while j<len(inequalities_list):
            k=j+1
            while k<len(inequalities_list):
                 if inequalities_list[j][i]*inequalities_list[k][i]==-1.0:
                     c=c+1
                 k=k+1
            j=j+1 
        count.append(c)    
        var.append(i)
        #print(count)
    #print(var)
        #print(count)
    #if len(count)!=0:  
    #print(count)
    min_index=count.index(min(count))        
    #min_index=min(enumerate(count), key=lambda x:x[1] if x[1]>0 else float('inf'))[0] #fastest I could get
    print(var[min_index], "is chosen variable")
    return var[min_index]   


def advantage2(inequalities_list, variables_to_eliminate):
    count=[]
    smallest_list=[]
    
    
    to_remove_unobserved=[]
    to_remove_unobserved.append(inequalities_list)
    args=[]
    for y in variables_to_eliminate:
           dummy=to_remove_unobserved[:]
           dummy.append(y)
           args.append(dummy)
    j=0 
    argument=[]     
    while j<len(args):
         argument.append(tuple(args[j]))
         j=j+1
         
    if __name__ == '__main__':
        with Pool(25) as pool:     
             smallest_list=list(pool.starmap(Fourier_Motzkin, argument))
    count=[]
    for i in smallest_list:
          count.append(len(i))
    min_index=count.index(min(count)) 
    print(variables_to_eliminate[min_index], "is the variable elimianted")
    print(len(smallest_list[min_index]),"is the length of non-redundant list after elimination of this variable")
    return (smallest_list[min_index],min_index);
    
    
    
    #to_remove_unobserved=inequalities_list[:]
    #for x in variables_to_eliminate:
         #unobserved_removed=[]
         #unobserved_removed.append(Fourier_Motzkin(to_remove_unobserved, x))
         #smallest_list.append(unobserved_removed[0])
         #count.append(len(unobserved_removed[0]))
    #min_index=count.index(min(count)) 
    #print(count)
    #print(min_index)
    #print(smallest_list[24])
    #return (smallest_list[min_index],min_index);
    
       
    
"To generate the inequalities corresponding to $ \mathcal{C} $"        
def C_cone_new(par,allnodes,obs,N):
    cond_indep=Markov_cond(par,obs,N)
    print(cond_indep, "for C cone")
    to_remove_unobserved=remove_redundancies(cond_indep,N,allnodes,obs)
   #print(to_remove_unobserved,len(to_remove_unobserved))
   #print("This is the C_cone after removing constraints....Looks good because its has some elements in decimals now.....and the remove constraint was atleast doing the calculations correctly" )
    #j=0
    #print(val2nice(allnodes & ~obs))
    to_remove_unobserved_copy=to_remove_unobserved[:]
    z,y,l=0,0,0
    while z<len(to_remove_unobserved):
        list_constraints=to_remove_unobserved_copy[:]
        list_constraints.pop(l)
        list_constraints=[[-1* i for i in j] for j in list_constraints]
        c=[0* i for i in range(len(list_constraints))]
        res=linprog_mosek(to_remove_unobserved[z], list_constraints, c)
        if res > -0.5:
             to_remove_unobserved_copy.remove(to_remove_unobserved[z])
        else:
            l=l+1
        z=z+1    
    to_remove_unobserved=to_remove_unobserved_copy[:]    
    print(to_remove_unobserved, len(to_remove_unobserved))
    print("This is the C_cone after removing redundant inequalities")  
    a=0
    variables_to_eliminate=[]
    res=[]
    for i in bitsof_inbits(allnodes & ~obs):
        j=i
        while j<= (2**N):
             if i in nonemptysubsetsof(j):
                 res.append(j)
                 a=a+1
             j=j+1
    b=0
    for y in res:
        if y not in variables_to_eliminate:
            variables_to_eliminate.append(y)
    length=len(variables_to_eliminate)   
    #print(variables_to_eliminate)
    while b<length:
        #print(variables_to_eliminate,"yoyoyoyo")
        (to_remove_unobserved,min_index)=advantage2(to_remove_unobserved,variables_to_eliminate)
        variables_to_eliminate.pop(min_index)
        b=b+1
    return to_remove_unobserved;       
        
    
 
    
"To generate the inequalities corresponding to $ \mathcal{C} $"    
def C_cone(par,allnodes,obs,N):
    cond_indep=Markov_cond(par,obs,N)
    print(cond_indep, "for C cone")
    to_remove_unobserved=remove_redundancies(cond_indep,N,allnodes,obs)
   #print(to_remove_unobserved,len(to_remove_unobserved))
   #print("This is the C_cone after removing constraints....Looks good because its has some elements in decimals now.....and the remove constraint was atleast doing the calculations correctly" )
    #j=0
    #print(val2nice(allnodes & ~obs))
    to_remove_unobserved_copy=to_remove_unobserved[:]
    z,y,l=0,0,0
    while z<len(to_remove_unobserved):
        list_constraints=to_remove_unobserved_copy[:]
        list_constraints.pop(l)
        list_constraints=[[-1* i for i in j] for j in list_constraints]
        c=[0* i for i in range(len(list_constraints))]
        res=linprog_mosek(to_remove_unobserved[z], list_constraints, c)
        if res > -0.5:
             to_remove_unobserved_copy.remove(to_remove_unobserved[z])
        else:
            l=l+1
        z=z+1    
    to_remove_unobserved=to_remove_unobserved_copy[:]    
    print(to_remove_unobserved, len(to_remove_unobserved))
    print("This is the C_cone after removing redundant inequalities")  
    a=0
    variables_to_eliminate=[]
    res=[]
    for i in bitsof_inbits(allnodes & ~obs):
        j=i
        while j<= (2**N):
             if i in nonemptysubsetsof(j):
                 res.append(j)
                 a=a+1
             j=j+1
    b=0
    for y in res:
        if y not in variables_to_eliminate:
            variables_to_eliminate.append(y)
    length=len(variables_to_eliminate)
    while b<length:
    #for i in bitsof_inbits(allnodes & ~obs):
          k=advantange1(to_remove_unobserved, variables_to_eliminate)
    #while j<= (2**N):
        #if i in nonemptysubsetsof(j):
          unobserved_removed=[]
          unobserved_removed.append(Fourier_Motzkin(to_remove_unobserved, k))
          to_remove_unobserved=unobserved_removed[0][:]
          variables_to_eliminate.remove(k)
          print("left variables", len(variables_to_eliminate))
          b=b+1    
     
    #print(to_remove_unobserved, len(to_remove_unobserved))
    #print("no")
    return to_remove_unobserved;  


def C_cone_check(par,allnodes,obs,N):
    cond_indep=Markov_cond(par,obs,N)
    print(cond_indep, "for C cone")
    to_remove_unobserved=remove_redundancies(cond_indep,N,allnodes,obs)
   #print(to_remove_unobserved,len(to_remove_unobserved))
   #print("This is the C_cone after removing constraints....Looks good because its has some elements in decimals now.....and the remove constraint was atleast doing the calculations correctly" )
    #j=0
    #print(val2nice(allnodes & ~obs))
    to_remove_unobserved_copy=to_remove_unobserved[:]
    z,y,l=0,0,0
    while z<len(to_remove_unobserved):
        list_constraints=to_remove_unobserved_copy[:]
        list_constraints.pop(l)
        list_constraints=[[-1* i for i in j] for j in list_constraints]
        c=[0* i for i in range(len(list_constraints))]
        res=linprog_mosek(to_remove_unobserved[z], list_constraints, c)
        if res > -0.5:
             to_remove_unobserved_copy.remove(to_remove_unobserved[z])
        else:
            l=l+1
        z=z+1    
    to_remove_unobserved=to_remove_unobserved_copy[:]    
    print(to_remove_unobserved, len(to_remove_unobserved))
    print("This is the C_cone after removing redundant inequalities") 
    a=0
    variables_to_eliminate=[]
    res=[]
    for i in bitsof_inbits(allnodes & ~obs):
        j=i
        while j<= (2**N):
             if i in nonemptysubsetsof(j):
                 res.append(j)
                 a=a+1
             j=j+1
    
    for y in res:
        if y not in variables_to_eliminate:
            variables_to_eliminate.append(y)
    length=len(variables_to_eliminate)
    #b=[24,28,20,22,16,18,17,21,40,42,3,14,0,1,2,4,5,6,7,34,8,9,10,11,12,13,15,19,41,23,43,25,26,27,29,30,31,32,45,33,35,36,38,37,39,44,46,47]
    #b=[49,57,41,45,33,37,35,43,34,42,7,29,1,3,5,9,11,13,15,10,17,19,21,23,25,27,31,39,38,47,46,51,53,55,59,61,63,2,54,6,14,18,26,22,30,50,58,62]
    #print(variables_to_eliminate[22])
    #b= [50, 51, 34, 35, 83, 3, 6, 127, 7, 22, 26, 46, 62, 82, 114, 84, 24, 40, 2, 10, 42, 43, 98, 115, 99, 11, 14, 15, 18, 19, 23, 27, 30, 31, 38, 39, 47, 54, 59, 63, 66, 67, 55, 58, 70, 71, 74, 75, 103, 102, 78, 79, 86, 87, 90, 91, 94, 95, 106, 107, 52, 20, 36, 41, 110, 111, 118, 116, 100, 119, 122, 123, 126, 4, 12, 5, 13, 21, 85, 28, 29, 44, 45, 60, 68, 69, 76, 92, 108, 109, 124, 8, 9, 56, 72, 88, 89, 25, 37, 53, 77, 101, 117, 73, 57, 61, 93, 125, 104, 105, 121, 120]
    #b=[50, 51, 34, 35, 83, 3, 6, 127, 7, 22, 26, 46, 62, 82, 114, 84, 24, 40, 2, 10, 42, 43, 98, 115, 99, 11, 14, 15, 18, 19, 23, 27, 30, 31, 38, 39, 47, 54, 59, 63, 66, 67, 55, 58, 70, 71, 74, 75, 103, 102, 78, 79, 86, 87, 90, 91, 94, 95, 106, 107, 52, 20, 36, 41, 110, 111, 118, 116, 100, 119, 122, 123, 126, 4, 12, 5, 13, 21, 28, 29, 44, 45, 37, 101, 53, 85, 117, 25, 60, 61, 68, 69, 76, 77, 92, 108, 109, 124, 8, 9, 56, 57, 72, 73, 88, 89, 93, 125, 104, 105, 120, 121]
    #b=[26, 90, 18, 82, 106, 3, 9, 17, 19, 83, 42, 58, 21, 85, 84, 20, 81, 28, 23, 25, 27, 33, 37, 43, 45, 87, 93, 101, 127, 22, 38, 50, 54, 52, 122, 12, 5, 7, 11, 13, 15, 31, 35, 39, 41, 105, 47, 29, 49, 113, 114, 51, 115, 92, 89, 53, 117, 116, 57, 121, 61, 1, 65, 67, 69, 71, 77, 91, 59, 97, 73, 75, 79, 55, 119, 63, 95, 99, 103, 109, 123, 107, 111, 125, 2, 6, 10, 14, 30, 34, 98, 46, 44, 62, 60, 66, 70, 74, 78, 86, 118, 94, 102, 110, 126, 4, 36, 68, 76, 100, 108, 124]
    for i in variables_to_eliminate: 
          print("eliminated variable is:", i)        
          unobserved_removed=[]
          unobserved_removed.append(Fourier_Motzkin(to_remove_unobserved, i))
          to_remove_unobserved=unobserved_removed[0][:]
          
    return to_remove_unobserved   
   
        

def perfectly_correlated_entropy_vector(allnodes,obs,N,correlated_nodes):
    particular_vector=[0 for i in range(2**N)]
    i=0
    while i<len(particular_vector):
        if ((i & correlated_nodes) in nonemptysubsetsof(correlated_nodes)) and (i & (allnodes & ~obs)==0):
                particular_vector[i]=1
        else:
                particular_vector[i]=0
        i=i+1        
    print(particular_vector)
    return particular_vector



def special_entropy_vector(allnodes,obs,N,random1,random2,xored,fixed):
    particular_vector=[0 for i in range(2**N)]
    i=0
    while i<len(particular_vector):
            if  (i & (allnodes & ~obs)==0):
                 if (random1|random2) in nonemptysubsetsof(i) or (random2|xored) in nonemptysubsetsof(i) or  (random1|xored) in nonemptysubsetsof(i):
                        particular_vector[i]=2
                        
                 else:  
                        particular_vector[i]=1
            else:
                particular_vector[i]=0
            i=i+1
    particular_vector[fixed]=0
    particular_vector[0]=0
    particular_vector[xored]=0
    particular_vector[fixed|xored]=0
    print(particular_vector)
    print(fixed)
    return particular_vector        




def interesting_special_entropy_vector(par,allnodes,obs,N,random1,random2,xored,fixed):
    random1=nice2val(random1)
    random2=nice2val(random2)
    xored=nice2val(xored)
    fixed=nice2val(fixed)
    vector=special_entropy_vector(allnodes,obs,N,random1,random2,xored,fixed)
    cond_indep=Markov_cond(par,obs,N)    
    to_remove_unobserved=remove_redundancies(cond_indep,N,allnodes,obs) 
    to_remove_unobserved=[[-1*i for i in j]for j in to_remove_unobserved]
    var_constraint,var_bound,i,component_sum=[],[0*j for j in range(len(to_remove_unobserved))],0,0
    while i<len(vector):
        var_constraint=[0*j for j in range(len(vector))]
        if i in nonemptysubsetsof(obs):
            var_constraint[i]=1
            to_remove_unobserved.append(var_constraint)
            var_bound.append(vector[i])
        component_sum= component_sum + vector[i]*vector[i]
        i=i+1
    vector=[-1*i for i in vector]
    res=-1*linprog_mosek(vector,to_remove_unobserved, var_bound)
    print(component_sum)
    if res<component_sum- 0.001:
       print(res) 
       print("This vector is interesting", vector) 
       return vector
    else:
        print(res)
        print("This vector is not interesting")
        return 0
    
    
                 
#interesting_special_entropy_vector([0b0,0b0,0b0,0b110,0b1001,0b11],0b111111,0b111100,6,'F', 'E','C','D')

    

def interesting_perfectly_correlated_entropy_vector(par,allnodes,obs,N,correlated_nodes):
    vector=perfectly_correlated_entropy_vector(allnodes,obs,N,correlated_nodes)
    cond_indep=Markov_cond(par,obs,N)    
    to_remove_unobserved=remove_redundancies(cond_indep,N,allnodes,obs) 
    to_remove_unobserved=[[-1*i for i in j]for j in to_remove_unobserved]
    var_constraint,var_bound,i,component_sum=[],[0*j for j in range(len(to_remove_unobserved))],0,0
    while i<len(vector):
        var_constraint=[0*j for j in range(len(vector))]
        if i in nonemptysubsetsof(obs):
            var_constraint[i]=1
            to_remove_unobserved.append(var_constraint)
            var_bound.append(vector[i])
        component_sum= component_sum + vector[i]*vector[i]
        i=i+1
    vector=[-1*i for i in vector]
    res=-1*linprog_mosek(vector,to_remove_unobserved, var_bound)
    print(component_sum)
    if res<component_sum- 0.001:
       print(res) 
       print("This vector is interesting", vector) 
       return vector
    else:
        print(res)
        print("This vector is not interesting")
        return 0    

#interesting_perfectly_correlated_entropy_vector([0b0,0b0,0b0,0b110,0b101,0b11],0b111111,0b111000,6,0b111000)           

#print([val2nice(i) for i in range(2**6)  ])
(['', '', '', 'BC', 'AC', 'ABE', 'ABDF'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'AD'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'ADF'], 'DEFG')
(['', '', '', 'C', 'BCD', 'ADE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'ABD'], 'DEFG')
(['', '', '', 'BC', 'AC', 'BE', 'ADF'], 'DEFG')
(['', '', '', 'BC', 'AD', 'ACE', 'ABF'], 'DEFG')
(['', '', '', 'AC', 'BC', 'AE', 'ABDF'], 'DEFG')
(['', '', '', 'AC', 'ABD', 'BC', 'A'], 'DEFG')
(['', '', '', 'C', 'BCD', 'ACE', 'AB'], 'DEFG')
(['', '', '', 'A', 'CD', 'BDE', 'BC'], 'AEFG')
(['', '', '', '', 'BCD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'BC', 'AD', 'CE', 'ABEF'], 'DEFG')
(['', '', '', 'BC', 'AD', 'CE', 'ABF'], 'DEFG')
(['', '', '', 'BC', 'C', 'ACDE', 'AB'], 'DEFG')
(['', '', '', 'BC', 'AD', 'ACE', 'ABEF'], 'DEFG')
(['', '', '', 'AC', 'BD', 'BCE', 'AEF'], 'DEFG')    
    

"""Fourier-Motzkin elimination eliminating variable at matrix[position] 
from the list of all inequalities in list matrix"""    
def Fourier_Motzkin(matrix, position):
    i, motzkin, j = 0, [], 0
    tech=0
    while i<len(matrix):
        j=i+1
        while j<len(matrix):
            #print("hello", matrix[j][position])
            if matrix[i][position] * matrix[j][position]<0.0 :#and round(abs(matrix[j][position]))!=0:
                #print("hello", matrix[j][position])
                motzkin.append(matrix_manipulate(matrix[i], matrix[j], position))
            j=j+1
        if matrix[i][position]==0:
            motzkin.append(matrix[i])    
        i=i+1 
    #motzkin=inequality_redundancy(motzkin)    
    motzkin1, k,l = motzkin[:], 0,1
    print(len(motzkin1), "length of redundant inequality")
    if tech==0:
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
              #print("linear program number=",l)
              start=time.time()
              #print(len(motzkin1))
              res = linprog_mosek(c, motzkin1, b)
              #print("time taken=", time.time()-start)
              l=l+1
              #print(motzkin1)
              #print(res['fun'], res['success'], res['status'], res['message'])
              if res > -0.5:
                  #print(motzkin[k])
                  motzkin.pop(k)
                  motzkin1=motzkin[:]
                  #print(res)
                  #print(res['fun'], res['success'])
              else:
                  #print(res,"kkkkkkkkkkkkkkk")
                  k=k+1
                  motzkin1=motzkin[:]
              for z in motzkin1:
                  if all(v==0.0 for v in z):
                      motzkin.remove(z)
                      motzkin1=motzkin[:]
          else:
              k=k+1;
      print(len(motzkin), "the lenght after removing redundancies")   
      return motzkin;
    else:
        motzkin2=inequality_redundancy(motzkin)
        return motzkin2

    
   
    
#Working - tests Fourier_Motzkin and matrix_manipulate functions
#a=[[0,0,1,5],[0,1,0,6],[1,0,0,4],[0,1,1,2]]
#c=0 #shuffle the above values according to the comment in green below to see the effect
#d=Fourier_Motzkin(a,c)
#print(d)
#print(len(d))
    
    

"""When the number of inequalities are quite large a VERY BIG ADVANTAGE CAN be
got from eliminating that variable in each step for which the  total number of pairs
of inequalities in which the variable occurs with opposite signs is the least.
another advantage is in eliminating thatfor index, m in enumerate(constraints_copy[0]):
                                               if m!=0 and index!=0 and (allnodes & ~obs) in nonemptysubsetsof(index):
                                                       equalities.append(matrix_manipub=0
while b<nonemptysubsetsof(15):
    print(b)
    b=b+1late(constraints_copy[j],constraints_copy[0], index))
                                                       break; variable for which the maximum number
of inequalities has its corresponding coefficient 0"""    


                     
#Testing Markov_cond- Working (with strat==add for all the dgas I tested)
#print(Markov_cond([0b0,0b1,0b10,0b110,0b1101,0b11101,0b110111],0b1011011,7))



def inequality_redundancy(motzkin):
    motzkin1,i=[],0
    while i<len(motzkin):
        if i==0:
             motzkin1.append(motzkin[i])
        else:
            c=motzkin[i] 
            b=[0 for i in range(len(motzkin1))]
            motzkin2=[[-1* i for i in j] for j in motzkin1]
            res = linprog_mosek(c, motzkin2, b)
            if res < -0.5:
                motzkin1.append(c)
        i=i+1        
    return motzkin1;
        
            
            
            
        


"To list inequalities in $\mathcal{C}$ but not in $mathcal\{I}"
def gdag_cone(par,allnodes,obs,N):
    list_I_cone=I_cone(par,obs,N) 
    print(list_I_cone, len(list_I_cone))
    print("This is the I_cone after removing constraints")
    list_C_cone=C_cone_check(par,allnodes,obs,N)
    print(list_C_cone, len(list_C_cone))
    print("This is the C_cone after Fourier Motzkin elimination of all unobserved variables' occurences")
    interesting_inequalities, non_interesting=[], []
    list_I_cone_reverse=[[-1*i for i in j]for j in list_I_cone]
    b=[0*j for j in range(len(list_I_cone))]
    for i in list_C_cone:
         res = linprog_mosek(i, list_I_cone_reverse , b)
         if res > -0.5:
             non_interesting.append(i)     
         else:
            interesting_inequalities.append(i)
   #print(non_interesting, len(non_interesting))
   #print("This is the list  of non interesting inequalities that are in C and I") 
    return  interesting_inequalities 


a=gdag_cone([0b0,0b0,0b11,0b100,0b1010,0b11001],0b111111,0b111100,6)  
#a=gdag_cone([0b0,0b0,0b11,0b100,0b1010,0b11001],0b111111,0b111100,6)  
#a=gdag_cone([0b0,0b0,0b11,0b1],0b1111,0b1110,4)            
#a=gdag_cone([0b0,0b0,0b0,0b110,0b1000,0b10101,0b110011],0b1111111,0b1111000,7)
(['', '', '', 'BC', 'AC', 'ABE', 'ABDF'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'AD'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'ADF'], 'DEFG')
(['', '', '', 'C', 'BCD', 'ADE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'DEFG')
(['', '', '', 'BC', 'AC', 'ABE', 'ABD'], 'DEFG')
(['', '', '', 'BC', 'AC', 'BE', 'ADF'], 'DEFG')
(['', '', '', 'BC', 'AD', 'ACE', 'ABF'], 'DEFG')
(['', '', '', 'AC', 'BC', 'AE', 'ABDF'], 'DEFG')
(['', '', '', 'AC', 'ABD', 'BC', 'A'], 'DEFG')
(['', '', '', 'C', 'BCD', 'ACE', 'AB'], 'DEFG')
(['', '', '', 'A', 'CD', 'BDE', 'BC'], 'AEFG')
(['', '', '', '', 'BCD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'C', 'BCD', 'ACDE', 'AB'], 'CEFG')
(['', '', '', 'BC', 'AD', 'CE', 'ABEF'], 'DEFG')
(['', '', '', 'BC', 'AD', 'CE', 'ABF'], 'DEFG')
(['', '', '', 'BC', 'C', 'ACDE', 'AB'], 'DEFG')
(['', '', '', 'BC', 'AD', 'ACE', 'ABEF'], 'DEFG')
(['', '', '', 'AC', 'BD', 'BCE', 'AEF'], 'DEFG')   

#[0b0,0b0,0b11,0b0,0b1010,0b10101],0b111111,0b111100,6) none of the 2 inequalities via different advantage same
#[0b0,0b0,0b0,0b1,0b1100,0b11010,0b110],0b1111111,0b1110001,7)
#[0b0,0b0,0b0,0b101,0b1011,0b110,0b1],0b1111111,0b1111000,7)
#[0b0,0b0,0b0,0b110,0b1001,0b10011],0b111111,0b111100,6) one inequality out of 2 via different advantage same
#[0b0,0b0,0b1,0b110,0b1100,0b10010],0b111111,0b111001,6) both inequalities via different advantage same
#[0b0,0b0,0b11,0b101],0b1111,0b1110,4) the inequality is same via both advantages
#[0b0,0b0,0b0,0b11,0b110],0b11111,0b11110,5) 0 inequalities via both advantages
#[0b0,0b0,0b11,0b110,0b101],0b11111,0b11100,5) the inequality is the same via both advantages
print(a,len(a),'hello') 
print("This is the list of inequalities that are in C_cone but not in I_cone") 
c=0  
        
#for i in a:
    #if i[1]==0 and i[3]==0 and i[5]==0 and i[7]==0 and i[9]==0 and i[11]==0 and i[13]==0 and i[15]==0:
        #print("coefficients of A, AB, AC, ABC, AD, ABD, ACD, ABCD all are 0")
        #c=c+1
   
    #else:
        #print("in", i, "coefficients of A, AB, AC, ABC, AD, ABD, ACD, ABCD all are not 0")
#print(c)
        



"To convert interesting inequalities from numbers to letters"
def interesting_inequalities(inequality,N):
    interesting_inequalities,b=[],2**N
    for i in inequality:
        c=[]
        j=0
        while j<b:
            c.append(str(i[j])+val2nice(j))
            j=j+1
        interesting_inequalities.append(c)    
        
    #print(a)
    return interesting_inequalities;  
print("List of Interesting Inequalities")      
print(interesting_inequalities(a,6))
print("--- %s seconds ---" % (time.time() - start_time))






               
