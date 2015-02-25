# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 10:33:05 2015

@author: Rafael
"""

def sdd(events,probs):
    import random
    nprobs=[x*1000 for x in probs] #so, here i multiply each float in 'probs' by 1000 and store the products in 'nprobs'
    newlist=[]
    for a in range(len(events)) : #then, in this loop, i create a list (newlist), in which each event appears 1000*its probability times
        b=nprobs[a]
        b=int(b)
        for c in range(b) :
            newlist.append(events[a]) 
    return (random.choice(newlist)) #and finally, i ramdonly sample an element from newlist
    

def contmark(states,q,v) :
    """IMPORTANT:rows and columns in the q-matrix MUST be in the order ACGT
    IMPORTANT2: the q-matrix should be a numpy matrix object, a 'list of lists' won't work
    IMPORTANT3: this function is dependendt on my discrete sampling function, 'sdd', pasted above"""
    import scipy as sp
    import random
    t=sp.linalg.expm(q*v) #get the transition probability matrix from the q-matrix
    chain=[]
    times=[]
    elapsedtime=0 #initialize a variable for the elapsed time
    currstate=random.choice(states) #randomly choose the initial state
    chain.append(currstate)
    wt=0 #initialize the waiting time variable
    while elapsedtime<v :
        """this series of if statements chooses the appropriate lambda (rate of the exponential distribution), given the currente state"""
        if currstate=='a' :
            lambd=-(q.item(0,0))
        if currstate=='c' :
            lambd=-(q.item(1,1))
        if currstate=='g' :
            lambd=-(q.item(2,2))
        if currstate=='t' :
            lambd=-(q.item(3,3))
        wt=random.expovariate(lambd) #use the appropriate lambda to draw a waiting time
        times.append(wt)
        elapsedtime+=wt
        """this series of if/elif statements draws the next state"""
        if currstate=='a' :
            currstate=sdd(states,t[0])
        elif currstate=='c':
            currstate=sdd(states,t[1])
        elif currstate=='g':
            currstate=sdd(states,t[2])
        elif currstate=='t':
            currstate=sdd(states,t[3])       
        chain.append(currstate)
    return chain, times

