# -*- coding: utf-8 -*-
"""
Created on Tue Apr 07 16:45:24 2015

@author: Rafael
"""

from EX6_CLEAN import Tree
from EX6_CLEAN import Node
from CTMCtree import ctmc 
import numpy
import scipy


          
def sitearray(node,siteindex):
    """takes as arguments a node and a site index, and outputs a prob array 
    (in order a,c,g,t) for that site on that node, based on child nodes, or, if node is a tip,
    on nucleotide present in its sequence"""
    """an attribute (a list) of markov object in each node is created, called sitearrays, that will contain prob
    arrays for each site in that node (thus it will be a list of lists)"""  
    if len(node.children)==0: #if node is a tip
        #print 'len(children) is 0, node is a tip'
        if node.seq[siteindex]=='a':
            #print 'nc at site at tip is a'
            array=[1,0,0,0]
        if node.seq[siteindex]=='c':
            #print 'nc at site at tip is c'
            array=[0,1,0,0]
        if node.seq[siteindex]=='g':
            #print 'nc at site at tip is g'
            array=[0,0,1,0]
        if node.seq[siteindex]=='t':
            #print 'nc at site at tip is t'
            array=[0,0,0,1]

        node.markov.sitearrays[siteindex]=array
        return array
    else: #if node is not a tip
        """print 'node is not a tip. calculations for child 0:' """
        #print 'starting calculatiosn for child 0'
        pmatrix=scipy.linalg.expm(node.children[0].markov.q*node.children[0].brl)#get pmatrix based on child's brl
        #print 'pmatrix for child 0', pmatrix
        
        probaa=pmatrix[0,0]#prob of a->a transition along branch linking current node to child 0
        probac=pmatrix[0,1]#prob of a->c transition
        probag=pmatrix[0,2]# etc
        probat=pmatrix[0,3]
        probachild0=(probaa*node.children[0].markov.sitearrays[siteindex][0]+ #multiplies transition probs by prob that child has respective nc at site, and sums up to get prob of a at this site at this node, considering only child 0
        probac*node.children[0].markov.sitearrays[siteindex][1]+
        probag*node.children[0].markov.sitearrays[siteindex][2]+
        probat*node.children[0].markov.sitearrays[siteindex][3]) 
        #print 'probachild0=', probachild0
    
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild0=(probca*node.children[0].markov.sitearrays[siteindex][0]+
        probcc*node.children[0].markov.sitearrays[siteindex][1]+
        probcg*node.children[0].markov.sitearrays[siteindex][2]+
        probct*node.children[0].markov.sitearrays[siteindex][3]) 
        #print 'probcchild0=',probcchild0
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild0=(probga*node.children[0].markov.sitearrays[siteindex][0]+
        probgc*node.children[0].markov.sitearrays[siteindex][1]+
        probgg*node.children[0].markov.sitearrays[siteindex][2]+
        probgt*node.children[0].markov.sitearrays[siteindex][3])
        #print 'probgchild0=',probgchild0
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild0=(probta*node.children[0].markov.sitearrays[siteindex][0]+
        probtc*node.children[0].markov.sitearrays[siteindex][1]+
        probtg*node.children[0].markov.sitearrays[siteindex][2]+
        probtt*node.children[0].markov.sitearrays[siteindex][3])  
        #print 'probtchild0=',probtchild0
                      
        """calculations for child 1:"""
        #print 'starting calculatiosn for child 1'
        pmatrix=scipy.linalg.expm(node.children[1].markov.q*node.children[1].brl)#get pmatrix based on child's brl
        #print 'pmatrix for child 0', pmatrix
        
        probaa=pmatrix[0,0]
        probac=pmatrix[0,1]
        probag=pmatrix[0,2]
        probat=pmatrix[0,3]
        probachild1=(probaa*node.children[1].markov.sitearrays[siteindex][0]+
        probac*node.children[1].markov.sitearrays[siteindex][1]+
        probag*node.children[1].markov.sitearrays[siteindex][2]+
        probat*node.children[1].markov.sitearrays[siteindex][3]) 
        #print 'probachild1=',probachild1
        
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild1=(probca*node.children[1].markov.sitearrays[siteindex][0]+
        probcc*node.children[1].markov.sitearrays[siteindex][1]+
        probcg*node.children[1].markov.sitearrays[siteindex][2]+
        probct*node.children[1].markov.sitearrays[siteindex][3]) 
        #print 'probcchild1=',probcchild1
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild1=(probga*node.children[1].markov.sitearrays[siteindex][0]+
        probgc*node.children[1].markov.sitearrays[siteindex][1]+
        probgg*node.children[1].markov.sitearrays[siteindex][2]+
        probgt*node.children[1].markov.sitearrays[siteindex][3])
        #print 'probgchild1=',probgchild1
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild1=(probta*node.children[1].markov.sitearrays[siteindex][0]+
        probtc*node.children[1].markov.sitearrays[siteindex][1]+
        probtg*node.children[1].markov.sitearrays[siteindex][2]+
        probtt*node.children[1].markov.sitearrays[siteindex][3])  
        #print 'probtchild1=',probtchild1

        """multiply probs for children 0 and 1 and get final array"""
       # print 'calculating array'
        array=[probachild0*probachild1,probcchild0*probcchild1,probgchild0*probgchild1,probtchild0*probtchild1]
        node.markov.sitearrays[siteindex]=array
        return array            

def findtips(node,siteindex,tiplist=None):
    """takes a node as argument and returns a list containing all descendent nodes that do not yet have
    sitearrays associated to them. In more technical terms, i suppose that is to say that it finds the tips of
    the current tree after it has been pruned of the subtrees whose site prob arrays have alraedy been calculated"""


    if tiplist == None:
        tiplist=[]
    if len(node.children)==0 and len(node.markov.sitearrays[siteindex])!=4: #if node is a (real) tip and does not yet have array for that site
        #print 'reached a tip:',node.name, 'Appending it to list'
        tiplist.append(node)
        #print 'list so far is:',tiplist
    
    elif len(node.markov.sitearrays[siteindex])!=4 and len(node.children[0].markov.sitearrays[siteindex])==4 and len(node.children[1].markov.sitearrays[siteindex])==4 : #second and third conditions assure that it is the last node without array for that site
        #print 'reached a last node without sitearrays',node.name,'Appending it to list'
        tiplist.append(node)
        #print 'list so far is:',tiplist
    else:
        for child in node.children:
            findtips(child,siteindex,tiplist=tiplist)
    return tiplist


def sitearraytree(node,siteindex):
    """takes a node (typically the root) root as argument and uses function sitearray (above) to calculate arrays for a given
    site at every node in tree."""
    """ the idea is to use my findtips function (above), that outputs a list of the last (closer to tips) nodes that do 
    not yet have site arrays, and then use my function sitearray to calculate arrays for nodes in that list. That is 
    implemented as a while loop that goes on until the root node itself has a site array for that site- that is, until
    node.sitearrays[siteindex] becomes a list with for itens (4 probs, corresponding to each nc)."""
    """IT WORKS!!!!!"""
    while len(node.markov.sitearrays[siteindex])!=4 : 
        for tip in findtips(node,siteindex) :
            sitearray(tip,siteindex)
      
        
def treelike(root):
    """uses sitearraytree (above) to calculate arrays for every site in seq, then sums itens in each array
    to get site likelihoods, and finally multiplies these to get tree likelihood"""
    for site in range(root.markov.nsites):
        sitearraytree(root,site)
    finalarrays=root.markov.sitearrays
    #multiply by the PIs:
    for array in finalarrays:
        array[0]=array[0]*root.markov.piarray[0]
        array[1]=array[1]*root.markov.piarray[1]
        array[2]=array[2]*root.markov.piarray[2]    
        array[3]=array[3]*root.markov.piarray[3]
    print finalarrays[0]
    print '---'
    sitelikes=[]
    for i in range(len(finalarrays)):
        sitelikes.append(sum(finalarrays[i]))
    print sitelikes[0]
    print '---'
    treelikelihood=1
    for item in sitelikes:
        treelikelihood*=item
    return treelikelihood
