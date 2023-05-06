#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 21:29:47 2022

@author: sk1808
"""
import time
import mosek
from linear_programming import linprog_mosek
from multiprocessing import Pool


start_time = time.time()

"""To choose the variable to eliminate at each step which generates the least redundant list of inequalities"""
def advantange1(inequalities_list, variables_to_eliminate):
    var=[]
    count=[]
    for i in variables_to_eliminate:
        j=0
        c=0
        while j<len(inequalities_list):
            k=j+1
            while k<len(inequalities_list):
                 if inequalities_list[j][i]*inequalities_list[k][i]<0.0:
                     c=c+1
                 k=k+1
            j=j+1 
        count.append(c)    
        var.append(i)
    min_index=count.index(min(count))        
    #min_index=min(enumerate(count), key=lambda x:x[1] if x[1]>0 else float('inf'))[0] #fastest I could get
    #print(var[min_index], "is chosen variable")
    return var[min_index]   



"""To eliminate the variable at each step which generates the least non-reduddant list of inequalities"""
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
        with Pool(125) as pool:     
             smallest_list=list(pool.starmap(Fourier_Motzkin, argument))
    count=[]
    for i in smallest_list:
          count.append(len(i))
    min_index=count.index(min(count)) 
    print(variables_to_eliminate[min_index], "is the variable elimianted")
    print((smallest_list[min_index]),"is the non-redundant list after elimination of this variable")
    return (smallest_list[min_index],min_index);



"""To round the elements of the matrix up-to desired accuracy"""
def matrix_round(matrix):
  round_matrix=[]  
  for i in matrix:
        if abs(i)<0.000000001:
            round_matrix.append(0)
        else:
            round_matrix.append(i)
  return round_matrix;  



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
    k=matrix1[position]/matrix2[position]
    if k>0.0: 
        matrix=[matrix1[j] - (abs(k) * matrix2[j]) for j in range(len(matrix1))]
        return matrix_round(matrix);
    if k<0.0:
        matrix=[matrix1[j] + (abs(k) * matrix2[j]) for j in range(len(matrix1))]                      
        return matrix_round(matrix);
    if k==0.0:
        return matrix_round(matrix1); 



"""Fourier-Motzkin elimination eliminating variable at matrix[position] 
from the list of all inequalities in list matrix"""    
def Fourier_Motzkin(matrix, position):
    i, motzkin, j = 0, [], 0
    #print(matrix,position)
    while i<len(matrix):
        j=i+1
        while j<len(matrix):
            if matrix[i][position] * matrix[j][position]<0.0 :
                motzkin.append(matrix_manipulate(matrix[i], matrix[j], position))
            j=j+1
        if matrix[i][position]==0:
            motzkin.append(matrix[i])    
        i=i+1    
    motzkin1, k,l = motzkin[:], 0,1
   #print(motzkin,"after FM list")
    while k<len(motzkin):
          start=0.0
          c=[i * 1 for i in motzkin1[k]]
          motzkin1.pop(k)
          if len(motzkin1)!=0:
              b=[0 for i in range(len(motzkin1))]
              motzkin1=[[-1* i for i in j] for j in motzkin1]
              #start=time.time()
              res = linprog_mosek(c, motzkin1, b)
              l=l+1
              if res > -0.5:
                  motzkin.pop(k)
                  motzkin1=motzkin[:]
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



"""Faster Fourier-Motzkin elimination using advantage1 techinque"""
def faster_elimination_1(inequalities_list, variables_to_eliminate):
    for i in range(len(variables_to_eliminate)):
          k=advantange1(inequalities_list, variables_to_eliminate)
          unobserved_removed=[]
          unobserved_removed.append(Fourier_Motzkin(inequalities_list, k))
          inequalities_list=unobserved_removed[0][:]
          variables_to_eliminate.remove(k)
    return (inequalities_list, len(inequalities_list));



"""Faster Fourier-Motzkin elimination using advantage2 techinque"""
def faster_elimination_2(inequalities_list, variables_to_eliminate):
    for i in range(len(variables_to_eliminate)):
         (inequalities_list, min_index)=advantage2(inequalities_list, variables_to_eliminate)
         variables_to_eliminate.pop(min_index)
    return inequalities_list;


print("--- %s seconds ---" % (time.time() - start_time))
"Testing-works fine"
inequalities=matrix=[[1,0,0,-1,-1,-1,0,0,0],[2,1,0,0,2,0,0,-1,0],[2,0,1,0,0,2,0,0,-1], [-1,0,-2,0,0,-3,0,0,1],[0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[3,-1,-1,2,0,0,-1,0,0],[-5,-1,-1,-2,-2,-2,1,1,1],[0,0,-1,0,0,-1,0,2,0]]
var=[1]
#print(faster_elimination_1(inequalities,var))