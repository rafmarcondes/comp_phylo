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
    

class MarkovChain(object): #define class
    def __init__(self, matrix,staspa): #define the variables i want associated with this instance of class
        self.matrix=matrix
        self.staspa=staspa
    def getchain(staspa,matrix,length) : #method to run a simulation and get a chain
        import random
        chain=[]
        currstate=random.choice(staspa)
        chain.append(currstate)
        for i in range ((length-1)) :
            for j in range (len(staspa)) :
                if currstate==staspa[j] :
                    currstate=sdd(staspa,matrix[j])
                    chain.append(currstate)
                    break
        return(chain)
        #end of class definition
    

mymarkov=MarkovChain(matrix=[[0.4,0.6],[0.50,0.5]],staspa=("A","B")) #create an instance of MarkovChain class  

print mymarkov.matrix #Try to retrieve 'matrix' variable from object - it works!

print mymarkov.getchain() #Try to run a method from the object - it is not working. Error message says i need two more arguments, but
#i thought that since the 'matrix' and 'staspa' arguments are variables associated with this object, i wouldnt
#need to pass them onto the method expicitly
