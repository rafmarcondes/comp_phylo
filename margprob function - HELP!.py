# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 14:43:02 2015

@author: Rafael
"""

def margprob(initial,final,q,v) :
    """A function to find the margprob of a final state, given the initial state. Rows and lines on q-matrix
    must be in order a, c, g, t."""
    import scipy    
    pmatrix=scipy.linalg.expm(q*v)
    iindex=0
    findex=0
    if initial=='a':
        iindex=0
    elif initial=='c':
        iindex=1
    elif initial=='g':
        iindex=2
    elif initial=='t':
        iindex=3
    if final=='a' :
        findex==0
    elif final=='c' :
        findex==1
    elif final=='g' :
        findex==2
    elif final=='t' :
        findex==3
    mp=pmatrix[iindex,findex]
    return mp


import numpy
qmatrix=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355')

print(margprob('t','a',qmatrix,0.1))

"""The function runs with no error message, but the output is not as expected. The q-matrix I am using is the one
 on Huelsenback page 45, and he gives corresponding p-matrices for several values of v on page 52, so I can 
 check if my function is working correctly. And it is not. The value it spits as output is a value from the p-matrix,
 but not from the cell it should be. This is driving me crazy, I've been trying to find out what is wrong with 
 this function since the afternoon yesterday. I'm sure I'm just missing something obvious, but can anyone 
 help me please?!"""



