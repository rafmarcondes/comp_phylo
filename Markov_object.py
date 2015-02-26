# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 10:01:11 2015

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
    

class MarkovChain(object): #define class. whatever is inside the parentheses tells python what class markov inherits 
#attributes from. (object) tells that it inherits from nothing
    def __init__(self, q,states, v): #a constructor. defines the variables i want associated with this 
    #instance of class
        self.q=q #a q-matrix
        self.states=states
        self.v=v
       
    def contmark(states,q,v) : #method to run a simulation and get a chain
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
        return chain

        #end of class definition
        

""""ok, created the class, now ket's try to use it:"""
import numpy
mymarkov=MarkovChain (q=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355'), states=('a','c','g','t'), v=5) #create an instance of MarkovChain class  

print mymarkov.q #Try to retrieve matrix variable from object - it works!

print mymarkov.contmark #NOT WORKING
