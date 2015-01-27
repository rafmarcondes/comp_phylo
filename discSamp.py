# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 16:03:48 2015

@author: Rafael

This is code written to complete a discrete sampling exercise for Jeremy Brown's CompPhylo class at LSU, Spring 2015.
Comments preceded by "INSTRUCTIONS" are Jeremy's literal instructions for the exercise; all other comments are mine.
"""

"""
INSTRUCTIONS: (1) Write a function that multiplies all consecutively decreasing numbers between a maximum and 
a minimum supplied as arguments. (Like a factorial, but not necessarily going all the way to 1).
"""

def cm(min,max) :
#the function is named after 'consecutive multiplication'. Its arguments are the minimum and maximum values. The
#output is the result of the multiplication, and it is not automatically printed to the screen
    result=max #since ranges in python are exclusive of the greater value, I initialize the product of the multiplication
    # with it, as it will not be multiplied in the for loop below.
    for i in range (min,max) :
        result=result*i
    return result
    
"""
INSTRUCTIONS: (2) Using the function you wrote in (1), write a function that calculates the binomial 
coefficient (see Definition 1.4.12 in the probability reading). Actually, do this twice. The first time (2a) 
calculate all factorials fully. 
"""

def bc(n,k) :
#this function takes as input two integers, n and k, in which k<n, and returns as output their binomial coefficient,
#that is, the number of possible subsets with size k, from a set with size n.
    coef=((cm(1,n))/((cm(1,(n-k)))*(cm(1,k))))
    return coef
    
"""
INSTRUCTIONS: Now re-write the function and cancel as many terms as possible so you can avoid unnecessary 
multiplication (see the middle expression in Theorem 1.4.13
"""


def fbc(n,k) :
#this function (named for Fast Binomial Coefficient) also takes as input two integers, n and k, in which k<n,
# and also returns as output their binomial coefficient

#this bit of code creates a list containing all integers n, n-1, n-2,...,n-k+1 (the factors in the
#numerator in the middle expression in Theorem 1.4.13)
    nlist=[]
    for i in range ((n-k+1),(n+1)):
        nlist.append(i)

#then, this bit of code multiplies all members of that list 
    p=1
    for j in nlist:
        p=p*j
    
#and finally, the result is divided by k!, to obtain the binomial coefficient
    bc=(p/(cm(1,k)))
    return bc
    
"""
INSTRUCTIONS: (3) Try calculating different binomial coefficients using both the functions 
from (2a) and (2b) for different values of n and k. Try some really big values there is a noticeable 
difference in speed between the (2a) and (2b) function. Which one is faster? By roughly how much?
"""

#print(bc(100000,50000)) #6.66sec
#print(fbc(100000,50000)) #3.78sec

#the 2b function is almost twice as fast

"""
INSTRUCTIONS: (4) Use either function (2a) or (2b) to write a function that calculates the probability of k 
successes in n Bernoulli trials with probability p. This is called the Binomial(n,p) distribution. 
See Theorem 3.3.5 for the necessary equation.
"""

def ber(k,n,p) :
    """this function calculates the probability of obtaining k sucesses in n Bernoulli trials (a trial in
which the only two possible outcomes are sucess or failure), in each of which the probability of success
is p. The output of the function is the compound probability and is named pp in the context of the function."""
    pp=(fbc(n,k))*(pow(p,k))*(pow((1-p),(n-k)))
    return pp
    
"""
INSTRUCTIONS:
(5) Now write a function to sample from an arbitrary discrete distribution. This function should take two 
arguments. The first is a list of arbitrarily labeled events and the second is a list of probabilities 
associated with these events. Obviously, these two lists should be the same length.
"""

def sdd(events,probs):
    """This function takes as arguments two lists, 'events' and 'probs', in which 'probs' contains floats 
    representing probabilities associated, respectively, with the itens in
    'events', that can be strings, floats, ints, or any data type. The two lists must be of same size.
    
    The strategy for this function is to randomly sample from a list (newlist) in which each event (elements in 'events') 
    appears a numebr of 
    times proportional to its probability (elements in 'probs')  """
    
    import random
    nprobs=[x*1000 for x in probs] #so, here i multiply each element in 'probs' by 1000 and store the products in 'nprobs'
    newlist=[]
    for a in range(len(events)) : #then, in this loop, i create a list (newlist), in which it event appears 1000*its probability times
        b=nprobs[a]
        b=int(b)
        for c in range(b) :
            newlist.append(events[a]) 
    return (random.choice(newlist)) #and finally, i ramdonly sample an element from newlist
    
"""INSTRUCTIONS: Imagine that you have a multiple sequence alignment with two kinds of sites. One type of site 
pattern supports the monophyly of taxon A and taxon B. The second type supports the monophyly of taxon A and taxon C.

(6) For an alignment of 400 sites, with 200 sites of type 1 and 200 of type 2, sample a new alignment 
(a new set of site pattern counts) with replacement from the original using your function from (5). Print out the 
counts of the two types.
"""

align=(['a']*200+['b']*200)
prob=[0.0025]*400 #0.0025 being 1/400, assuming equal probabilities of sampling each site

newalign=[]
for n in range (400):
    newalign.append((sdd(align,prob)))
    
acount=0
bcount=0
for n in newalign :
    if n=='a' :
        acount+=1
    if n=='b':
        bcount+=1
        
print('In the new alignment, there were', acount, 'a-sites, and', bcount, 'b-sites')

"""INSTRUCTIONS: (7) Repeat (6) 100 times and store the results in a list."""

align=(['a']*200+['b']*200)
prob=[0.0025]*400 #0.0025 being 1/400, assuming equal probabilities of sampling each site

countapertrial=[]
countbpertrial=[]
for m in range(100):
    newalign=[]
    for n in range (400):
        newalign.append((sdd(align,prob)))
    
    acount=0
    bcount=0
    for n in newalign :
        if n=='a' :
            acount+=1
        if n=='b':
            bcount+=1
    countapertrial.append(acount)
    countbpertrial.append(bcount)
    
print(countapertrial)
   

"""INSTRUCTIONS: (8) Of those 100 trials, summarize how often you saw particular proportions of type 1 vs. type 2."""

"""I use a histogram to summarize the number of trials in which I obtained particular counts of a-sites.
I do not plot a histogram for counts of b-sides because that would be redundant, since b=400-a"""

import matplotlib.pyplot as plt
plt.hist(countapertrial, bins=100)
plt.xlabel('number of a-sites')
plt.ylabel('number of trials')

"""INSTRUCTIONS: (9) Calculate the probabilities of the proportions you saw in (8) using the binomial 
probability mass function (PMF) from (4)."""

"""i'm not sure i understand the instructions correctly. i assume that the goal is to 
calculate the probability of obtaining the proportion of a-sites (or b-sites) that i saw in (6) - not (8)"""

"""i saw 211 a-sites (successes), out of 400 trials. therefore:"""
print(ber(211,400,0.5))

"""INSTRUCTIONS: (10) Compare your results from (8) and (9).

(11) Repeat 7-10, but use 10,000 trials."""


