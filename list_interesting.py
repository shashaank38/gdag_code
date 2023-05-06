#!/usr/bin/env python2
"abcdefgh"
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags

# Settings
MAXSIZE = 5
NUMPROCESSES = 7 # number of parallel processes

# List interesting DAGs of size at most MAXSIZE
#for n in range(1, MAXSIZE+1):
a=enumdags.enumdags(6, NUMPROCESSES)

for x in a:
    y=x
    print(y[len(y)-1])
    print(x)



#for x in enumdags.subsetsof(obs):
        #for y in enumdags.subsetsof(obs & ~x):
             #for z in enumdags.subsetsof(obs & ~x & ~y):
                  #w= obs & ~x & ~y & ~z
                  #if(x|y|z|w==obs):
                       #if(eseparation_functions.theorem(x,y,z,w, par,obs)==True):
                           #print([[bitfuncs.val2nice(i) for i in par], obs])
                           #print("Is interesting")
                           #break
             #else:
                 #continue
             #break
        #else:
            #continue
        #break
    #print(([[bitfuncs.val2nice(i) for i in par], obs]))
    #print("Could not be checked")