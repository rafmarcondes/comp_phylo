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
    def __init__(self, q, states,v): #a constructor. defines the variables i want associated with this 
    #instance of class
        self.q=q #a q matrix
        self.states=states #the state space
        self.v=v #a branch length
      
        
    def contmark(self) :        
        """IMPORTANT:rows and columns in the q-matrix MUST be in the order a,c,g,t
        IMPORTANT2: the q-matrix should be a numpy matrix object, a 'list of lists' won't work
        IMPORTANT3: this function is dependendt on my discrete sampling function, 'sdd', pasted above"""
        import scipy as sp
        import random
        statup=tuple(self.states) 
        chain=[]
        times=[]
        elapsedtime=0 #initialize a variable for the elapsed time
        currstate=random.choice(self.states) #randomly choose the initial state
        chain.append(currstate)
        wt=0 #initialize the waiting time variable
        while elapsedtime<self.v :
            if currstate=='a' :
                lambd=-(self.q.item(0,0))
            if currstate=='c' :
                lambd=-(self.q.item(1,1))
            if currstate=='g' :
                lambd=-(self.q.item(2,2))
            if currstate=='t' :
                lambd=-(self.q.item(3,3))
            wt=random.expovariate(lambd) #use the appropriate lambda to draw a waiting time
            times.append(wt)
            elapsedtime+=wt
            if currstate=='a' :
                self.states.remove('a') #in order not to draw the same state as the current one
                currstate=sdd(self.states,(self.q.item(1),self.q.item(2),self.q.item(3)))
            elif currstate=='c' :
                self.states.remove('c')
                currstate=sdd(self.states,(self.q.item(4),self.q.item(6),self.q.item(7)))
            elif currstate=='g' :
                self.states.remove('g')
                currstate=sdd(self.states,(self.q.item(8),self.q.item(9),self.q.item(11)))
            elif currstate=='t' :
                self.states.remove('t')
                currstate=sdd(self.states,(self.q.item(12),self.q.item(13),self.q.item(14))   )
            chain.append(currstate)
            self.states=list(statup) #reset the 'states' list to its original composition
        return chain, times
    

""""ok, created the class, now let's try to use it:"""

import numpy
bla=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355')
mymarkov=MarkovChain(q=bla, states=['a','c','g','t'], v=5) #create an instance of MarkovChain class  

print mymarkov.q #retrieve matrix variable from object
print mymarkov.contmark() #use contmark method to obtain a chain
