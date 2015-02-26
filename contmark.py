# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 10:33:05 2015

@author: Rafael
"""

def sdd(events,probs):
    """This function is named after "Sample from Discrete Distribution". It takes as arguments two lists,
    'events' and 'probs', in which 'probs' contains floats 
    representing probabilities associated, respectively, with the itens in 'events', that can be strings, 
    floats, ints, or any data type. The two lists must be of same size."""
    
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
    statup=tuple(states)
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
            states.remove('a')
            currstate=sdd(states,(q.item(1),q.item(2),q.item(3)))
        elif currstate=='c' :
            states.remove('c')
            currstate=sdd(states,(q.item(4),q.item(6),q.item(7)))
        elif currstate=='g' :
            states.remove('g')
            currstate=sdd(states,(q.item(8),q.item(9),q.item(11)))
        elif currstate=='t' :
            states.remove('t')
            currstate=sdd(states,(q.item(12),q.item(13),q.item(14))   )
        chain.append(currstate)
        states=list(statup)
    return chain, times
    
