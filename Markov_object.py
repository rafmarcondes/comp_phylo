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
    def __init__(self, q, states,v,chain,times): #a constructor. defines the variables i want associated with this 
    #instance of class
        self.q=q 
        self.states=states
        self.v=v
        self.chain=chain
        self.times=times
      
            
    def contmark(self) :
        """IMPORTANT:rows and columns in the q-matrix MUST be in the order ACGT
        IMPORTANT2: the q-matrix should be a numpy matrix object, a 'list of lists' won't work
        IMPORTANT3: this function is dependendt on my discrete sampling function, 'sdd', pasted above"""
        import scipy as sp
        import random
        statup=tuple(self.states) 
        """during the function, i'll have to remove elements from the 'states' list, but later reset the list 
        to its original composition. i'll use the 'statup' tuple as a 'backup'to do that"""
        chain=[]
        times=[]
        elapsedtime=0 #initialize a variable for the elapsed time
        currstate=random.choice(self.states) #randomly choose the initial state
        chain.append(currstate)
        wt=0 #initialize the waiting time variable
        while elapsedtime<self.v :
            """this series of if statements chooses the appropriate lambda (rate of the 
            exponential distribution), given the currente state"""
            if currstate=='a' :
                lambd=-(self.q.item(0,0))
            if currstate=='c' :
                lambd=-(self.q.item(1,1))
            if currstate=='g' :
                lambd=-(self.q.item(2,2))
            if currstate=='t' :
                lambd=-(self.q.item(3,3))
            wt=random.expovariate(lambd) #use the appropriate lambda to draw a waiting time
            elapsedtime+=wt
            if elapsedtime>self.v : 
                """if elapsedtime will exceed v with the addition of the last wt, 
                then undo its addition to elapsed time and break the while loop, because v has been reached"""
                elapsedtime-=wt
                break
            times.append(wt)
            """this series of if/elif statements draws the next state"""       
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
        times.append(self.v-elapsedtime) 
        """append to the times list the difference between the elapsed time so far andv, so that the final elapsedtime, or sum(times), equals v"""
        self.chain=chain #updtate the attributes chain and times to contain the output of the simulation
        self.times=times
        return chain, times
        
    def timeperstate(self) :
        """OUTPUT: a list with 4 itens, each being the total time spent in states a, c, g and t, in that order"""
        ttime=0
        atime=0
        ctime=0
        gtime=0
        for n in range(0,len(self.chain)):
                if self.chain[n]=='a' :
                    atime+=self.times[n]
                if self.chain[n]=='c' :
                    ctime+=self.times[n]
                if self.chain[n]=='t' :
                    ttime+=self.times[n]
                if self.chain[n]=='g' :
                    gtime+=self.times[n]
        actgtimes=[]
        actgtimes.append(atime)
        actgtimes.append(ctime)
        actgtimes.append(gtime)
        actgtimes.append(ttime)
        return actgtimes
    
           


""""ok, created the class, now let's try to use it:"""


import numpy
bla=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355')

mymarkov=MarkovChain(q=bla, states=['a','c','g','t'], v=5, chain=[], times=[]) #create an instance of MarkovChain class  

print mymarkov.q #Try to retrieve matrix variable from object - it works!

print mymarkov.contmark() #this also wroks!!!! 

print mymarkov.timeperstate() #same here!!
