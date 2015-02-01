# -*- coding: utf-8 -*-
"""
An Introduction to Likelihood

@author: jembrown
"""

"""
There are two primary ways to use probability models. Given what we know to be true about an 
experimental setup, we can make predictions about what we expect to see in an upcoming trial. 
For this purpose, probability functions are what we need. If R is the outcome of an experiment (i.e., an event) 
and p is a parameter of the probability function defined for these experimental outcomes, we can write these 
expectations or predictions as:

P(R|p). - RSM: prediction of R (a result or outcome) given p (a parameter or probability). p is fixed and we
 want to know about R



This conditional probability tells us what to expect about the outcomes of our experiment, given knowledge 
of the underlying probability model. Thus, it is a probability function of R given a value for p (and the 
model itself).

However, we might also wish to ask what we can learn about p itself, given outcomes of trials that have 
already been observed. This is the purview of the likelihood. Likelihoods are functions of 
parameters (or hypotheses, more generally) given some observations. The likelihood function of a 
parameter value is defined as:

L(p;R) = P(R|p)

RSM: L(p;R) - likelihood of p, given R (a result, an observation, an outcome)

Note that this is the same probability statement we saw above. However, in this context we are 
considering the outcome (R) to be fixed and we're interested in learning about p. Note that the
 likelihood is sometimes written in several different ways: L(p;R) or L(p) or L(p|R). P(R|p) gives
 a probability when R is discrete or a probability density when R is continuous. Since likelihoods 
 are only compared for some particular R, we do not need to worry about this distinction. Technically 
 speaking, likelihoods are just said to be proportional to P(R|p), with the constant of proportionality
 being arbitrary.

There are some very important distinctions between likelihoods and probabilities. First, likelihoods do 
NOT sum (or integrate) to 1 over all possible values of p. Therefore, the area under a likelihood curve 
is not meaningful, as it is for probability.

It does not make sense to compare likelihoods across different R. For instance, smaller numbers of observations 
generally produce higher values of P(R|p), because there are fewer total outcomes.

Likelihood curves provide useful information about different possible values of p. When we are interested
 in comparing discrete hypotheses instead of continuous parameters, the likelihood ratio is often used:

L(H1;R)     P(R|H1)
-------  =  -------
L(H2;R)     P(R|H2)

Now, let's try using likelihoods to learn about unknown aspects of the process that's producing some data.


---> Inferring p for a binomial distribution <---

First, we'll start by trying to figure out the unknown probability of success associated with
 a Binom(3,p) random variable. If you want to try this on your own later, the following code will perform 
 draws from a binomial with 5 trials. You can simply change the associated value of p to whatever you'd like. 
 To make the inference blind, have a friend set this value and perform the draws from the Binomial for you,
 without revealing the value of p that they used.
"""

from scipy.stats import binom

n = 5
p = 0.5 # Change this and repeat

data = binom.rvs(n,p)

"""
For the in-class version of this exercise, I'm going to perform a manual draw from a binomial using colored
 marbles in a cup. We'll arbitrarily define dark marbles as successes and light marbles as failures.

Record the outcomes here:

Draw 1:
Draw 2:
Draw 3:
Draw 4:
Draw 5:

all marbles drawn were dark

Number of 'successes': 4

Now record the observed number of succeses as in the data variable below.
"""

data =4   # Supply observed number of successes here.
numTrials = 5

"""
Since we are trying to learn about p, we define the likelihood function as;

L(p;data) = P(data|p)

If data is a binomially distributed random variable [data ~ Binom(5,p)]

P(data=k|p) = (5 choose k) * p^k * (1-p)^(n-k)

So, we need a function to calculate the binomial PMF. Luckily, you should have just written one
 and posted it to GitHub for your last exercise. Copy and paste your binomial PMF code below. For now,
 I will refer to this function as binomPMF(). 
"""

#from JMB_binomial import * # I've stored my function in the file JMB_binomial.py
#import scipy

def cm(min,max) :
#the function is named after 'consecutive multiplication'. Its arguments are the minimum and maximum values. The
#output is the result of the multiplication, and it is not automatically printed to the screen
    result=max #since ranges in python are exclusive of the greater value, I initialize the product of the multiplication
    # with it, as it will not be multiplied in the for loop below.
    for i in range (min,max) :
        result=result*i
    return result
    
def fbc(n,k) :
    nlist=[]
    for i in range ((n-k+1),(n+1)):
        nlist.append(i)
    p=1
    for j in nlist:
        p=p*j

    bc=(p/(cm(1,k)))
    return bc
        
def ber(k,n,p) :
    """this function calculates the probability of obtaining k sucesses in n Bernoulli trials (a trial in
which the only two possible outcomes are sucess or failure), in each of which the probability of success
is p. The output of the function is the compound probability and is named pp in the context of the function."""
    pp=(fbc(n,k))*(pow(p,k))*(pow((1-p),(n-k)))
    return pp
    
    
    
"""
Now we need to calculate likelihoods for a series of different values for p to compare likelihoods. There 
are an infinite number of possible values for p, so let's confine ourselves to steps of 0.05 between 0 and 1.

""""

# Set up a list with all relevant values of p

pvalues=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]

# Calculate the likelihood scores for these values of p, in light of the data you've collected


ls=[]
for x in pvalues:
    p=ber(4,5,x)
    ls.append(p)

print (ls)
    


# Find the maximum likelihood value of p (at least, the max in this set)


ind=ls.index(max(ls))
print("the value of p with the max likelihood is", pvalues[ind],"its likelihood being", max(ls))

import matplotlib.pyplot as plt
plt.scatter(pvalues,ls)
plt.xlabel('p')
plt.ylabel('likelihood score')


# What is the strength of evidence against the most extreme values of p (0 and 1)?

#Total. their likelihood is zero.


# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)

lratio=[]
for l in ls :
    lr=l/(max(ls))
    lratio.append(lr)
print (lratio)

import matplotlib.pyplot as plt
plt.scatter(pvalues,lratio)
plt.xlabel('p')
plt.ylabel('likelihood ratio')


#MY GUESS: THE TRUE VALUE OF p IS IN THE INTERVAL (0.75,0.85)


"""
Now let's try this all again, but with more data. This time, we'll use 20 draws from our cup of marbles.
"""

data =12   # Supply observed number of successes here.
numTrials = 20


# Calculate the likelihood scores for these values of p, in light of the data you've collected

ls=[]
for x in pvalues:
    p=ber(12,20,x)
    ls.append(p)

print (ls)

# Find the maximum likelihood value of p (at least, the max in this set)

ind=ls.index(max(ls))
print("the value of p with the max likelihood is", pvalues[ind],"its likelihood being", max(ls))


import matplotlib.pyplot as plt
plt.scatter(pvalues,ls)
plt.xlabel('p')
plt.ylabel('likelihood score')


# What is the strength of evidence against the most extreme values of p (0 and 1)?

#Total. Their likelihoods are zero.


# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)

lratio=[]
for l in ls :
    lr=l/(max(ls))
    lratio.append(lr)
print (lratio)

import matplotlib.pyplot as plt
plt.scatter(pvalues,lratio)
plt.xlabel('p')
plt.ylabel('likelihood ratio')

# When is the ratio small enough to reject some values of p?
