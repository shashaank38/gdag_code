"""Main code for classifying GDAGS"""

# see Joe Henson, Raymond Lal, Matthew F. Pusey
# "Theory-independent limits on correlations from generalised Bayesian
#  networks", arXiv:1405.2572

# Copyright (C) 2014 Matthew F. Pusey -- see LICENSE.txt

import itertools, multiprocessing
from bitfuncs import *
from closure import *


cdef bint trystrat(int n, int obs, strat, substrat, dec_nodes, np,
        int not_dec_nodes, downstream, iobs):
    """Try a specific application of Theorem 26, namely steps 2-3 in
    Appendix C

    n -- number of nodes
    obs -- bitstring specifying which nodes are observable
    strat -- choices of root notes, in same order as dec_nodes
    substrat -- ordering for strat and dec_nodes we want to attempt
    dec_nodes -- "tricky nodes", in same order as strat
    np -- bitstring of parents for each node with step 1 already applied
    not_dec_nodes -- bitstring of observable nodes not in dec_nodes
    downstream -- for each node X, bitstring of observable nodes Y with
        X squiggly arrow to Y, i.e. path from X to Y through unobserved
    iobs -- observable conditional independences for this GDAG

    In the notation of Appendix C:
        The node number of T_i is dec_nodes[substrat[i]]
        The node number of R_i is strat[substrat[i]]

    Returns True if application has succeded and so C=I for this GDAG
    """
    cdef int targ, targbit, calc, consider, x, xind, ipar, ind, i, ibit
    s_np = np[:]
    cdef int already = 0

    # Step 2
    for i in substrat:
        targ = dec_nodes[i]
        targbit = 1<<targ
        calc = strat[i]

        s_np[targ] &= not_dec_nodes | already # Step 2(a)
        for j in bitsof(downstream[calc] & ~already & ~targbit):
            if s_np[targ] & ~s_np[j] == 0:
                s_np[j] |= targbit # Step 2(b)
        already |= targbit

    # [ Step 3 is implicit - we simply ignore the unoberved nodes ]

    # Calculate the standard generating set of conditional independences
    # (node indep. non-descendants given parents for each node) and see
    # if they are all observable independences for the original GDAG
    nondesc = [obs & ~(1<<x) for x in range(n)]
    for i in bitsof(obs):
        ibit = 1<<i
        consider = s_np[i]
        while consider:
            x = consider & ~(consider - 1)
            consider &= ~x
            xind = getindex(x)
            consider |= s_np[xind]
            nondesc[xind] &= ~ibit

    for i in bitsof(obs):
        ipar = s_np[i]
        ind = nondesc[i] & ~ipar
        if ind and Triplet(1<<i, ind, ipar) not in iobs:
            return False

    return True

def isclassical(int n, par, int obs):
    """ Determine if our sufficient condition tells us C = I

    (n, par, obs) is the standard format for a GDAG:
    n -- number of nodes
    par[i] -- bitstring of parents of node i
    obs -- bitstring of observable nodes
    """

    idag = dag_indeps(par)
    cdag = closure(idag)
    iobs = observable_indeps(cdag, obs)

    cdef int i, ibit, u_roots, ipar, i_is_obs, consider
    cdef int x, xind, xpar, xpar_u
    
    # Identify "tricky nodes" and the possible associated "root nodes"
    np = list(par)
    dec_nodes = []
    dec_options = []
    knows = []
    downstream = [0 for x in xrange(n)]
    not_dec_nodes = 0
    for i in xrange(n):
        ibit = 1<<i
        u_roots = 0
        ipar = par[i]
        i_is_obs = ibit & obs

        # Explore unobserved parents, and their unobserved parents etc
        consider = ipar & ~obs
        while consider:
            x = consider & ~(consider - 1)
            xind = getindex(x)
            xpar = par[xind]
            xpar_u = xpar & ~obs
            consider &= ~x
            consider |= xpar_u

            if i_is_obs:
                downstream[xind] |= ibit
            if xpar_u == 0:
                u_roots |= x
            # Step 1 of Appendix C can be applied once and for all
            np[i] |= xpar & obs

        if i_is_obs:
            if u_roots:
                dec_nodes.append(i)
                dec_options.append(bitsof(u_roots))
            else:
                not_dec_nodes |= ibit

    # If all observed nodes have only observed parents, trivially C=I
    if dec_nodes == []:
        return True

    # Now the two things we need to search over:
    # "every possible ordering of T"
    substrats = list(itertools.permutations(xrange(len(dec_nodes))))
    # "each element T_i associated with every possible R_i"
    for strat in itertools.product(*dec_options):
        for substrat in substrats:
            if trystrat(n, obs, strat, substrat, dec_nodes, np,
                    not_dec_nodes, downstream, iobs):
                return True

    return False

def implies(implicants, implicand): # implicand -> d, implicants-> trips
    """Does (c.i Triple) implicand follow from one of the implicants?"""
    for have in implicants:
        if implicand < have:
            return True
    return False

def observable_indeps(gen, int obs):
    """Calculate all the conditional independences that follow from gen   
    and only mention nodes in the bitstring obs"""
    res = []
    add=[]
    cdef int i, j, k
    for i in nonemptysubsetsof(obs):
        for j in nonemptysubsetsof(obs & ~i):
            for k in subsetsof(obs & ~i & ~j):
                ijk = Triplet(i,j,k)
                if implies(gen, ijk):
                    res.append(ijk)
                    add.append(i)
                    add.append(j)
                    add.append(k)
                    
                    
    return res #it was "return tuple(res)" here

def dag_indeps(par):
    """Calculate a generating set of conditional independces for a DAG

    par[i] -- bit array giving parents of i, assumed to all be before i
    """
    strat = []
    for i in xrange(len(par)):
        C = par[i]
        if (1<<i)-1 & ~C:
            strat.append(Triplet(1<<i, (1<<i)-1 & ~C, C))
    return strat;
    

cdef inline int applyperm(perm, int x):
    """Apply the permutation perm to the bitstring x"""
    cdef int y = 0, i, j
    for i,j in enumerate(perm):
        if ((1<<j) & x):
            y |= (1<<i)
    return y

cdef inline bint donethis(int n, par, int obs, done):
    """Is a GDAG isomorphic to (n, par, obs) in done?"""
    cdef int new, old, oldbit, newbit, newobs, newnode, oldnode
    newpartemplate = [0 for newnode in xrange(n)]
    for i in itertools.permutations(xrange(n)):
        newpar = newpartemplate[:]
        newobs = 0
        for new,old in enumerate(i):
            oldbit = 1<<old
            newbit = 1<<new
            for newnode, oldnode in enumerate(i):
                if oldbit & par[oldnode]:
                    newpar[newnode] |= newbit
            if oldbit & obs:
                newobs |= newbit

        if (tuple(newpar), newobs) in done:
            return True
    return False

def irreducible(n, par, obs):
    """Is (n, par, obs) worth listing after considernig reducibility?
    
    Returns False if is (n, par, obs) reducible to a smaller GDAG for
    which either
    (a) C and I are unchanged (i.e. reduction not by 1, 4, and 5 of
        Appendix D.1), or
    (b) Our sufficient condition fails for the smaller GDAG, and so the
        smaller GDAG will already be on the list of "interesting" ones
    """

    # Appendix D.1, 4: If there is a proof that C != I without some
    # obvervable node, then making that node trivial gives a proof for
    # par
    for i in bitsof(obs):
        noti = ~(1<<i)
        thispar = [x & noti for x in par]
        thispar[i] = 0
        if not isclassical(n, thispar, obs): # case (b) above
            return False

    # Appendix D.1, 5: Causal conditionals for nodes whose parents are
    # all observable can be seen in statistics. Hence we can enforce the 
    # uselessness" of an incoming edge and then use a proof without that
    # edge
    all_p_obs = 0
    for i in bitsof(obs):
        if par[i] & ~obs == 0:
            all_p_obs |= (1<<i)

    thispar = list(par)
    for i in bitsof(all_p_obs):
        for j in propersubsetsof(par[i]):
            thispar[i] = j
            if not isclassical(n, thispar, obs): # case (b) above
                return False
        thispar[i] = par[i]

    # Appendix D.2, 1: "Heisenberg picture" / composition of channels:
    # if a node has an unobserved parent with no other children, then
    # the channel can be incorporated into the measurement / composed
    # with channel
    for i in xrange(n):
        for j in bitsof(par[i] & ~obs):
            jbit = 1<<j
            if any([(jbit & parents and node != i) for 
                    node, parents in enumerate(par)]):
                continue
            return False # case (a) above

    # Appendix D.2, 2: If we have a parentless unobserved node with only
    # two children, one of which is observed and has no other parents
    # then all that unobserved node can achieve is deciding that child
    # might just use that child directly instead, giving smaller DAG
    u_nodes = ((1<<n) - 1) & ~obs
    for i in bitsof(u_nodes):
        if par[i]:
            continue
        ibit = 1<<i
        children = [j for j in xrange(n) if (par[j] & ibit)]
        if len(children) != 2:
            continue
        for j in children:
            if (1<<j) & obs and par[j] == ibit:
                return False # case (a) above

    # Appendix D.1, 6: Unobserved nodes whos parents and children are
    # subsets of another can be subsumed into that
    for i in bitsof(u_nodes):
        ibit = 1<<i
        for j in bitsof(u_nodes):
            if i == j:
                continue

            if par[i] & ~par[j]:
                continue

            jbit = 1<<j
            kidsok = True
            for k in xrange(n):
                if k != j and par[k] & ibit and not par[k] & jbit:
                    kidsok = False
                    continue

            if kidsok:
                return False # case (a) above

    return True

def trydag(args):
    """Return args=(n, par, obs) if an irreducibily interesting GDAG"""
    if not isclassical(*args) and irreducible(*args):
        return args

def trydag_count(args):
    """Does args=(n, par, obs) give an interesting GDAG?"""
    return isclassical(*args)

def find_candidates(int n, bint irrcheck):
    """Returns a list of GDAGs with n nodes that are worth looking at

    This means that:
    (a) No GDAG in the list is isomorphic to another
    (b) If irrcheck is True, all GDAGs pass initial reducibilty checks:
        (i)   They are connected
        (ii)  They don't have childless unobserved nodes
        (iii) They don't have unobserved nodes with a sole parent,
              also unobserved

    It is worth checking for a few reducibilities now since it avoids
    expensive isomorphism tests for obviously uninteresting GDAGs"""
    cdef int numbits = n*(n-1) // 2
    cdef int everything = (1<<n) - 1
    done = set()
    cdef int maxi = 1<<numbits
    cdef int i, j, ignore, mask, childless = 0
    cdef int connectedA, oldconnectedA, ipar, obs
    trythese = set()

    for i in xrange(maxi):
        par = []
        for j in xrange(n):
            ignore = j*(j-1) // 2
            mask = (1<<j) - 1
            par.append((i>>ignore) & mask)
        par = tuple(par)

        if irrcheck:
            childless = everything
            for i in par:
                childless &= ~i

            # Appendix D.1, 1: disconnected DAGs don't tell us more than
            # each component does
            connectedA = 1
            oldconnectedA = 0
            while connectedA != oldconnectedA:
                oldconnectedA = connectedA
                for i, ipar in enumerate(par):
                    if ipar & connectedA:
                        connectedA |= ipar | (1<<i)
            if connectedA != everything:
                continue

        for obs in xrange(1<<n):
            if irrcheck:
                # Appendix D.1, 2: Childless unobserved nodes are
                # pointless
                if childless & ~obs:
                    continue

                # Appendix D.1, 3: "Schrodinger picture": if an
                # unobserved node only has one parent, and its
                # unobserved, just apply channel to parent
                if any([par[i] & (par[i] - 1) == 0 and par[i] & ~obs for
                        i in bitsof(everything & ~obs)]):
                    continue

            if donethis(n, par, obs, done):
                continue
            done.add((par,obs))

            trythese.add((n, par, obs))

    return trythese

def enumdags(int n, int numprocesses):
    """Main interface to this module: output a list of interesting
    GDAGs of size n, using numprocesses processes in parallel"""
    trythese = find_candidates(n, True)
    
    res = map(trydag, trythese)
    for x in res:
        if x is not None:
            (n, par, obs) = x
            interesting=[[x for x in par],obs] 
            #print(interesting)
            yield interesting
            #print([val2nice(x) for x in par], val2nice(obs))
            
def countdags(int n, int numprocesses):
    """Alternative interface to this module: output n, the number of
    GDAGs, and the number of interesting ones, using numprocesses
    processes in parallel"""
    trythese = find_candidates(n, False)
    
    print(n, len(trythese))

"Try enumdags and countdags"
#enumdags(4,1)
#countdags(4,2)

""" Partially Wrong(earlier integer inputs were taken in decimal rather than 
binary"""
par=[0b0,0b1,0b10,0b110]

strat=dag_indeps(par)
#print(strat)
n=4
obs=0b1110
res=observable_indeps(strat, obs)
#print(res)
args=[4,par, 1]
a=trydag(args)
print(a)

"Check ALL the important functions for 5 Nodes Bell dag"
#n=5
#par=[0b0,0b0,0b0,0b11, 0b110]

#strat=dag_indeps(par)
#print(strat)
#print(" ")
#print(" ")
#obs=0b11101
#trips= closure(strat)
#print(trips)
#print(" ")
#print(" ")
#res=observable_indeps(trips, obs)
#print(res)
#print(" ")
#print(" ")

#b=isclassical(n, par,  obs)
#print(b)

"""Check ALL the important functions for 4 nodes Instrumental dag for which 
C=I criteria does not apply to reach any conclusion"""

#n=4
#par=[0b0,0b0,0b11,0b101]
#strat=dag_indeps(par)
#print(strat)
#print(" ")
#print(" ")
#obs=0b1110
#trips= closure(strat)
#print(trips)
#print(" ")
#print(" ")
#res=observable_indeps(trips, obs)
#print(res)
#print(" ")
#print(" ")

#b=isclassical(n, par,  obs)
#print(b)

"""Check ALL the important functions for 6 nodes Instrumental dag for which 
C=I criteria does not apply to reach any conclusion"""
"""n=7
par=[0b0,0b0,0b0,0b110, 0b101, 0b10010, 0b101001]

strat=dag_indeps(par)
print(strat)
print(" ")
print(" ")
obs=0b1111000
trips= closure(strat)
print(trips)#trips is not a list while strat and res are lists
print(" ")
print(" ")
res=observable_indeps(trips, obs)
print(res)
#print(" ")
#print(" ")"""

#b=isclassical(n, par,  obs)
#print(b)
#n=4
#par=[0b0,0b0,0b11,0b101]



#par.pop(2)
#par.insert(2,None)
#strat=dag_indeps(par)
#print(par)

"""n=3
par=[0b0, 0b1, 0b0]
obs=0b111
strat=dag_indeps(par)
#print(strat)
#print(" ")
trips=closure(strat)
#print(trips)
res=observable_indeps(trips, obs) 
for i in res:
  print("The observed conditional independencies in the original gdag are: ")
  print(i)
print(res)"""
#print(len(res)) 

"""j=0b10000
v=4#The index of the node to be deleted
k=val2nice(j)
print("node to be deleted is")
print(k)"""

"""Detect all the possible places where the deleted node comes any where 
in the string at that index"""
"""
def subgraph(par, int obs):
    parents=[]
    modified=[]
    i=0
    while i<len(par): #for x in par
      parents.append(val2nice(par[i])) #x instead of par[i]  
      i=i+1
    print("The original dag is")  
    print(parents)
    print("With observed nodes")
    print(val2nice(obs))
    i=0
    
    while i<len(parents):
        string1=parents[i]
        string2=''
        l=0
        while l<len(string1):
           if(string1[l]!=k):
              string2=string2+string1[l]
           if(string1[l]==k):
             string2=string2+''   
           l=l+1
        modified.append(string2)
        i=i+1
    modified[v]=''
    print("The subgraph is:") 
    print(modified)             
    new_par=[] 
    i=0
    while i<len(modified):
       new_par.append(nice2val(modified[i]))
       i=i+1
    print(new_par)
    return new_par

#new_par=subgraph(par,obs)

#print(new_par)
#strat=dag_indeps(par)
#print(strat)
#trips=closure(strat)
#print(trips)
#print(" ")
#res=observable_indeps(trips, obs)  
#print(res)   
a=input("Enter 1st observed set")
b=input("Enter 2nd observed set")

c=input("Enter the conditioning set")"""

"""
par=[0b0,0b0,0b11,0b0,0b1110,0b10101]#Jack dag
#par=[0b0,0b1,0b10,0b100,0b1000,0b10000]#chains
#par=[0b0,0b1,0b0,0b1010]#colliders
#par=[0b0,0b1,0b1,0b110,0b1001,0b1001]#collider fork

k='FB'
c='C'
def non_descendant():
    strat=dag_indeps(par)
                                                                                                                        trips=closure(strat)
    w=k
    z=c
    i=0
    n=[]
    summ=0
    while i<len(z):
        string1=z[i]
        d=nice2val(string1)
        j=0
        m=0
        while j<len(w):
            string2=w[j]
            e=nice2val(string2)
            position=ord(string2)-65
            string3=val2nice(par[position])
            print(string3)
            f=nice2val(string3)
            g=Triplet(d,e,f)
            if(string1 not in string3):
                if(implies(trips,g)==True):
                     m=m+0
                if(implies(trips,g)==False):    
                     m=m+1
                     break
            if(string1 in string3):
                  m=m+0
                 
            j=j+1
            print(m)
        n.append(m)
        i=i+1
    print(n)
    print(m)
    for x in n:
        summ=summ+x        
    print(summ)       
    if(summ==0):
       print("No node in Z is descendant from W")
    else:
       print("Some node in Z is descendant from W")         
            
            
non_descendant()            
"""

""

"""
def eseparation():
    i=nice2val(a)
    j=nice2val(b)
    k=nice2val(c)
    d=Triplet(i,j,k) # ask - if I put these 3 above conditions before the functions then output is wrong....why?? I have even changed the variables i,j,k to n,o,p...still the error
    new_par=subgraph(par,obs)
    print("To test e separartion in:")
    print(d)
    print(new_par)
    print(" ")
    strat=dag_indeps(new_par)
    print(strat)
    print(" ")
    trips=closure(strat)
    print(trips)
    print(" ")
    if(implies(trips,d)==True):
        print("True")
    else:
        print("False")
    #res=observable_indeps(trips, obs)  
    #print(res)
    #print(" ")
    #l=0
    #m=0
    #print(len(res))    
    #while l<len(res):
        #if(d==res[l]):
           #m=m+1
        #l=l+1
    #if(m!=0):
      #print("True")
    #else:
      #print("False") 
       
eseparation()  
#non_descendent()"""





            

    
          
        
    
    
    
    
    