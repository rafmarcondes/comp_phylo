# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 10:26:34 2015

@author: Rafael
"""

"""a function to find a maximum likelihood estimate of v"""

def margprob(initial,final,q,v) :
    """A function to find the margprob of a final state, given the initial state. Rows and lines on q-matrix
    must be in order a, c, g, t."""
    import scipy    
    pmatrix=scipy.linalg.expm(q*v)
    iindex=0
    findex=0
    if initial=='a':
        iindex=0
    if initial=='c':
        iindex=1
    if initial=='g':
        iindex=2
    if initial=='t':
        iindex=3
    if final=='a':
        findex=0
    if final=='c':
        findex=1
    if final=='g':
        findex=2
    if final=='t':
        findex=3
    mp=pmatrix[iindex,findex]
    return mp     


def mlv(initial,final,q,vStart,diff): 
    """a function to find the ML estimate of the branch length (v) between known initial and final states of a 
    continuous-time Markov chain. argument q is a q-matrix in the format if a numpy matrix. rows and columns must
    be in order a,c,g,t. vStart is initial value of v to be tried, and diff is initial step to next v to be tried"""
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
    
    vCurr=vStart
    
    while diff>0.001 :
        likeCurr=margprob(initial,final,q,vCurr)
        vUp=vCurr+diff
        vDown=vCurr-diff
        likeUp=margprob(initial,final,q,vUp)
        likeDown=margprob(initial,final,q,vDown)
        if likeDown>likeCurr :
            vCurr=vDown
            vUp=vCurr+diff
            vDown=vCurr-diff
        if likeUp>likeCurr :
            vCurr=vUp
            vUp=vCurr+diff
            vDown=vCurr-diff
        else :
            diff=(diff/2)            
    return vCurr
