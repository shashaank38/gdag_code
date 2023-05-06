#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 05:28:11 2022

@author: shashaank
"""

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
new_current_varaibles=[]    
def Fourier_Motzkin(matrix, position, histories, explicit_variables, implicit_variable):
    new_histories, new_explicit_variables, new_implicit_variables =[],[],[]
    new_current_varaibles.append(position)
    i, motzkin, j = 0, [], 0
    tech=0
    while i<len(matrix):
        j=i+1
        while j<len(matrix):
            if matrix[i][position] * matrix[j][position]<0.0 :#and round(abs(matrix[j][position]))!=0:
                matrix_manipulate_result=matrix_manipulate(matrix[i], matrix[j], position)
                motzkin.append(matrix_manipulate_result)
                new_histories.append([histories[i],histories[j]])
                new_explicit_variables.append( [explicit_variables[i], [position]])
                q=0
                implicit_var=[]
                while q<len(matrix_manipulate_result):
                    if q not in new_current_varaibles and matrix_manipulate_result[q]==0:
                        implicit_var.append(q)
                    q=q+1    
                new_implicit_variables.append(implicit_variables[i] + implicit_var)
            j=j+1   
        if matrix[i][position]==0:
            motzkin.append(matrix[i])  
            new_histories.append([histories[i]])
            new_explicit_variables.append(explicit_variables[i])
            new_implicit_variables.append([])
            
        #new_histories=list(dict.fromkeys(new_histories))               
        i=i+1 
    #motzkin=inequality_redundancy(motzkin)
    motzkin1, k,l = motzkin[:], 0,1   
    print(new_histories, "new histories")
    print(new_explicit_variables, "new explicit variables" )
    print(new_implicit_variables, "new implicit variables")
    print(new_current_varaibles)
    print(len(new_histories),len(new_explicit_variables), len(new_implicit_variables))
    #print(new_implicit_variables[2])
    print(motzkin, len(motzkin))
    #print(len(new_histories[4]), len(new_explicit_variables[4]),len(new_implicit_variables[4]))
    return 0

matrix=[[1,0,0,-1,-1,-1],[2,1,0,0,2,-1],[2,0,1,0,0,2], [-1,0,-2,0,0,-3],[0,0,1,0,0,0],[0,0,0,1,0,0],[3,-1,-1,2,0,0],[-5,-1,-1,-2,-2,-2],[0,0,-1,0,0,-1]]
position=1
histories=[i for i in matrix]
#print(histories)
explicit_variables=[[] for i in matrix]
#print(explicit_variables)
implicit_variables=[[] for i in matrix]
#print(implicit_variables)
current_variables=[]
print("after elimination of 1st variable")
Fourier_Motzkin(matrix, position, histories, explicit_variables, implicit_variables)
matrix=[[1, 0, 0, -1, -1, -1], [5.0, 0, -1.0, 2.0, 2.0, -1.0], [-3.0, 0, -1.0, -2.0, 0, -3.0], [2, 0, 1, 0, 0, 2], [-1, 0, -2, 0, 0, -3], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, -1, 0, 0, -1]]
histories=[[[1, 0, 0, -1, -1, -1]], [[2, 1, 0, 0, 2, -1], [3, -1, -1, 2, 0, 0]], [[2, 1, 0, 0, 2, -1], [-5, -1, -1, -2, -2, -2]], [[2, 0, 1, 0, 0, 2]], [[-1, 0, -2, 0, 0, -3]], [[0, 0, 1, 0, 0, 0]], [[0, 0, 0, 1, 0, 0]], [[0, 0, -1, 0, 0, -1]]]
explicit_variables=[[], [[], [1]], [[], [1]], [], [], [], [], []]
implicit_variables=[[], [[], []], [[], [4]], [], [], [], [], []]
position=2
print("after elimination of 2nd variable")
Fourier_Motzkin(matrix, position, histories, explicit_variables, implicit_variables)

   