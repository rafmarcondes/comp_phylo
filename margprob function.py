# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 14:43:02 2015

@author: Rafael
"""

def margprob(initial,final,q,v) :
    """A function to find the margprob of a final state, given the initial state. Q-matrix must be a 
    numpy matrix, and rows and lines must be in order a, c, g, t."""
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
        



