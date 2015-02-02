# -*- coding: utf-8 -*-
"""
Created on Sun Feb 01 15:52:49 2015

@author: Rafael
"""

"""
Sometimes it will not be feasible or efficient to calculate the likelihoods for every
value of a parameter in which we're interested. Also, that approach can lead to large
gaps between relevant values of the parameter. Instead, we'd like to have a 'hill
climbing' function that starts with some arbitrary value of the parameter and finds
values with progressively better likelihood scores. This is an ML optimization
function. There has been a lot of work on the best way to do this. We're going to try
a fairly simple approach that should still work pretty well, as long as our likelihood 
surface is unimodal (has just one peak). Our algorithm will be:

(1) Calculate the likelihood for our starting parameter value (we'll call this pCurr)
"""

def cm(min,max) :

    result=max 
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
    pp=(fbc(n,k))*(pow(p,k))*(pow((1-p),(n-k)))
    return pp
    
    
    
pCurr=0.5 # I came up with an arbitrary starting parameter value of 0.5 Its likelihood i calculated below:

likeCurr=ber(12,20,pCurr)
print(likeCurr)

"""
(2) Calculate likelihoods for the two parameter values above (pUp) and below (pDown)
our current value by some amount (diff). So, pUp=pCurr+diff and pDown=pCurr-diff. To
start, set diff=0.1, although it would be nice to allow this initial value to be set
as an argument of our optimization function.
"""
diff=0.1
pUp=pCurr+diff
pDown=pCurr-diff

likeUp=ber(12,20,pUp)
likeDown=ber(12,20,pDown)


"""
(3) If either pUp or pDown has a better likelihood than pCurr, change pCurr to this
value. Then repeat (1)-(3) until pCurr has a higher likelihood than both pUp and
pDown.
"""

if likeDown>likeCurr :
    pCurr=pDown
    pUp=pCurr+diff
    pDown=pCurr-diff
if likeUp>likeCurr :
    pCurr=pUp
    pUp=pCurr+diff
    pDown=pCurr-diff

"""
(4) Once L(pCurr) > L(pUp) and L(pCurr) > L(pDown), reduce diff by 1/2. Then repeat
(1)-(3).
"""

if likeDown>likeCurr :
    pCurr=pDown
    pUp=pCurr+diff
    pDown=pCurr-diff
if likeUp>likeCurr :
    pCurr=pUp
    pUp=pCurr+diff
    pDown=pCurr-diff
else :
    diff=(diff/2)

"""
(5) Repeat (1)-(4) until diff is less than some threshold (say, 0.001).
"""

likeCurr=ber(12,20,pCurr)

pUp=pCurr+diff
pDown=pCurr-diff

likeUp=ber(12,20,pUp)
likeDown=ber(12,20,pDown)

if likeDown>likeCurr :
    pCurr=pDown
    pUp=pCurr+diff
    pDown=pCurr-diff
if likeUp>likeCurr :
    pCurr=pUp
    pUp=pCurr+diff
    pDown=pCurr-diff
else :
    diff=(diff/2)

"""
(6) Return the final optimized parameter value.

I just ran the above code several times, until diff equaled 0.0004. The optimized parameter value then was 0.6

Write a function that takes some starting p value and observed data (k,n) for a
binomial as its arguments and returns the ML value for p.

To write this function, you will probably want to use while loops. The structure of
these loops is

while (someCondition):
    code line 1 inside loop
    code line 2 inside loop
    
As long as the condition remains True, the loop will continue executing. If the
condition isn't met (someCondition=False) when the loop is first encountered, the 
code inside will never execute.

If you understand recursion, you can use it to save some lines in this code, but it's
not necessary to create a working function.
"""

# Write a function that finds the ML value of p for a binomial, given k and n.


def findp(k,n,pStart,diff): #"""this is a function that uses a hill-climbing algorithm to find the parameter p (the probability of sucess) for a for a binomial, given k and n. The starting parameter and the initial step for the next parameter to be tried are provided by the user as arguments"""
    def cm(min,max) :
        result=max 
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
        pp=(fbc(n,k))*(pow(p,k))*(pow((1-p),(n-k)))
        return pp 
    
    pCurr=pStart
    
    while diff>0.001 :
        likeCurr=ber(k,n,pCurr)
        pUp=pCurr+diff
        pDown=pCurr-diff
        likeUp=ber(k,n,pUp)
        likeDown=ber(k,n,pDown)
        if likeDown>likeCurr :
            pCurr=pDown
            pUp=pCurr+diff
            pDown=pCurr-diff
        if likeUp>likeCurr :
            pCurr=pUp
            pUp=pCurr+diff
            pDown=pCurr-diff
        else :
            diff=(diff/2)            
    return pCurr
    

a=0
a=(findp(49,50,0.3,0.1))


"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of 
p. Now, we will empirically determine one way to construct such an interval. To do 
so, we will ask how far away from the true value of a parameter the ML estimate 
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then 
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once 
you have this distribution, find the likelihood ratio cutoff you need to ensure 
that the probability of seeing an LR score that big or greater is <= 5%. 
"""

# Set a starting, true value for p

trueP = 

# Simulate 1,000 datasets of 200 trials from a binomial with this p
# If you haven't already done so, you'll want to import the binom class from scipy:
# from scipy.stats import binom
# binom.rvs(n,p) will then produce a draw from the corresponding binomial.



# Now find ML parameter estimates for each of these trials



# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum
# likelihood (ML) in the denominator. Sort the results and find the value
# corresponding to the 95th percentile.



# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values.
# Find the 95th percentile of these values. Compare these values to this table:
# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look
# at the 0.05 column. Do any of these values seem similar to the one you calculated?
# Any idea why that particular cell would be meaningful?



# Based on your results (and the values in the table), what LR statistic value 
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?



# Using this cutoff, what interval might you report for the 5- and 20-trial data
# sets above?



# We've talked in previous classes about two ways to interpret probabilities. Which
# interpretation are we using here to define these intervals?



