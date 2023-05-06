"""Helper functions for dealing with bitstrings"""

# Copyright (C) 2014 Matthew F. Pusey
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.

import string

def getindex(int bit):
    """Get the position of the lowest 1 in bit"""
    cdef int i
    for i in xrange(26):
        if bit & 1<<i:
            return i

def val2nice(int val):
    """Output a bitstring val as letters, e.g. 101 becomes 'AC'"""
    return "".join([string.ascii_uppercase[i] for i in bitsof(val)])

def nice2val(nice):
    """Turn letters nice into a bitstring, e.g. 'AC' becomes 101"""
    cdef int val = 0
    for i in nice:
        val |= 1<<string.ascii_uppercase.find(i)
    return val

def bitsof(int left):
    """Generator of the locations of 1s in bitstring left"""
    cdef int x
    while left:
        x = left & ~(left - 1)
        left &= ~x
        yield getindex(x)
        
def bitsof_inbits(int left):
    """Generator of the bitstrings with a single 1 in the place of each 1 in bitstring left"""
    cdef int x
    while left:
        x = left & ~(left - 1)
        left &= ~x
        yield x        

def nonemptysubsetsof(int space):
    """Generator of the non-empty subsets of bitsring space"""
    cdef int x = space
    while x:
        yield x
        x = (x-1) & space
    #print(x)

def subsetsof(int space):
    """Generator of the subsets of bitsring space"""
    cdef int x
    for x in nonemptysubsetsof(space):
        yield x
    yield 0

def propersubsetsof(int space):
    """Generator of the proper subsets of bitsring space"""
    cdef int x = space
    while x:
        x = (x-1) & space
        yield x
    return x    
#val=0b100  
#b=val2nice(val)
#print(b)

#c=nice2val(b)
#print(bin(c))
#c='ACD'
#d=nice2val(c)
#e=bin(d)
#f=[e,e,e,e,e]
#print(f)

#W=0b10101
#v=16
#print(int(v))
#print(val2nice(W))
#print(val2nice(v))
#V=111011
#print(val2nice(W))
#for S  in subsetsof(W):
    #if v==S:
       #print("True")
#for S in subsetsof(W):    
    #print(val2nice(S))
#print(bin(59))
#print(bin(58))    

#print(nice2val("CAB"))

#w=16
#print(bin(w))
#print(val2nice(w))
#n=15
#print(bin(n))
#print(val2nice(n))
#n=n&~w
#print(bin(n))
#print(val2nice(n))
#i=0
#if i=1:
    #print("hey")
#else:
    #print("no")    
#z=0b1011
#i=2
#print(z & (1 << i))
#z=z&(1<<i)
#print(bin(z)) 
#z=z&~(1<<i)
#print(bin(z))

#for i in bitsof(z):
    #print(i)
    #for j in subsetsof(i):
        #print(val2nice(j))
#import math
#def countBits(number):
     
     #log function in base 2 
     #take only integer part
     #return int((math.log(number) /
                #math.log(2)) + 1);
 
# Driver Code
#num = 65;
#print(countBits(num));   

#a=0b000
#b=0b100
#if b in subsetsof(a):
    #print("y")
#else:
    #print("n") 
    """hi
print(val2nice(11),"hi",val2nice(7))    
for i in subsetsof(11):  
   print(val2nice(i))   
print(val2nice(11&7)) 
a=[]
for i in subsetsof(11):
    a.append(i)
print(a)
print(nice2val("D"))
a=[]
for i in range(0,16):
    a.append(val2nice(i))
print(a)    
    

  
 
print((1|2|4)) 

s=31
i=4
for j in (val2nice(s& ~i)):
    print(nice2val(j))
i=15
for j in subsetsof(i):
    print(val2nice(j))    
i=2
while i>2:
   print(i)
   i=i+1
for i in bitsof_inbits(31& ~8&~2):
   print(i) 
   
print(val2nice(31), val2nice(15))
print(31 & ~(1 << (4 - 1)))  
print(31&  ~16) 
print((1<<4)-1)
a=[1,2,3]
b=[4,5,6,7,8]
c=a+b
print(c)
a=[2,3,5]
b=[2,3,5]
c=[1,1,1]
print(a)
if(a==b):
    print("Y")
else:
    print("n") 
a.append(b)
a.append(c)
print(a)  
print(val2nice(14))
for i in bitsof_inbits(25):
    print(i)
print(bin(3))    

for i in bitsof_inbits(14):
    print(i)
for i in nonemptysubsetsof(14):
    print(val2nice)
a=[1,2,3,4,5]
b=[1,2,3]

print(a)

for i in range(1,1):
    print(i)
a=[[1,2,3],[2,3],[4,5,6,7],[8]]
a.remove([1,2])
print(a)
a,i=[1,2,3,4,5],0
while i<len(a):
    a.remove(a[i])
    print(a) 
    print(a[i])
    i=i+1
a=[[1,2,3],[4,5,6]]        
print(a[1][0])         
    
a=[1,2,3,4,5,6]
b=a[:]
for i in b:            
    a.remove(i)
    print(a)
a=[1,2,3,4,5,6]    
while a:
   i=a[0]
   a.remove(i)
   print(a)
   
a=[1,2,3]   
i=0
a.remove(3)
print(a)
a=[0 for i in range(10)]
print(a)
print(nice2val("ADE"))
a=31
b=25
print(~b)
print(a&~b)
print(bin(a&~b))
for i in bitsof_inbits(a&~b):
    print(i)
W=2
print(W & ~(1<<3)) 
i=7
if i in nonemptysubsetsof(3):
    print(i)
print(val2nice(14))
bound=(0,None)
bounds=[bound for i in range(10)]
print(bounds)
a=0
N=7
variables_to_eliminate=[]
for i in bitsof_inbits(7):
        j=i
        while j<= (2**N):
             if i in nonemptysubsetsof(j):
                 variables_to_eliminate.append(j)
                 a=a+1
             j=j+1
b=[]
variables_to_eliminate=list(set(variables_to_eliminate))
for j in variables_to_eliminate:
    b.append((j))
print(b) 
print(len(b)) ;
a=[]
print(min(a))"""  
print(nice2val('CE'))
