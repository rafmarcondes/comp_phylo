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

"""
class likelihood:"""
    

    
def sitearray(node,siteindex):
    """takes as arguments a node and a site index, and outputs a prob array 
    (in order a,c,g,t) for that site on that node, based on child nodes, or, if node is a tip,
    on nucleotide present in its sequence"""
    """an attribute (a list) of node object is created, called sitearrays, that will contain prob
    arrays for each site in that node (thus it will be a list of lists)"""  
    """this function seems to be working perfectly"""
    if len(node.children)==0: #if node is a tip
        #print 'len(children) is 0, node is a tip'
        if node.seq[siteindex]=='a':
            #print 'nc at site at tip is a'
            array=(1,0,0,0)
        if node.seq[siteindex]=='c':
            #print 'nc at site at tip is c'
            array=(0,1,0,0)
        if node.seq[siteindex]=='g':
            #print 'nc at site at tip is g'
            array=(0,0,1,0)
        if node.seq[siteindex]=='t':
            #print 'nc at site at tip is t'
            array=(0,0,0,1)
        node.sitearrays=[]
        node.sitearrays.append(array)
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
        probachild0=(probaa*node.children[0].sitearrays[siteindex][0]+ #multiplies tranition probs by prob that child has respective nc at site, and sums up to get prob of a at this site at this node, considering only child 0
        probac*node.children[0].sitearrays[siteindex][1]+
        probag*node.children[0].sitearrays[siteindex][2]+
        probat*node.children[0].sitearrays[siteindex][3]) 
        #print 'probachild0=', probachild0
    
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild0=(probca*node.children[0].sitearrays[siteindex][0]+
        probcc*node.children[0].sitearrays[siteindex][1]+
        probcg*node.children[0].sitearrays[siteindex][2]+
        probct*node.children[0].sitearrays[siteindex][3]) 
        #print 'probcchild0=',probcchild0
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild0=(probga*node.children[0].sitearrays[siteindex][0]+
        probgc*node.children[0].sitearrays[siteindex][1]+
        probgg*node.children[0].sitearrays[siteindex][2]+
        probgt*node.children[0].sitearrays[siteindex][3])
        #print 'probgchild0=',probgchild0
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild0=(probta*node.children[0].sitearrays[siteindex][0]+
        probtc*node.children[0].sitearrays[siteindex][1]+
        probtg*node.children[0].sitearrays[siteindex][2]+
        probtt*node.children[0].sitearrays[siteindex][3])  
        #print 'probtchild0=',probtchild0
                      
        """calculations for child 1:"""
        #print 'starting calculatiosn for child 1'
        pmatrix=scipy.linalg.expm(node.children[1].markov.q*node.children[1].brl)#get pmatrix based on child's brl
        #print 'pmatrix for child 0', pmatrix
        
        probaa=pmatrix[0,0]
        probac=pmatrix[0,1]
        probag=pmatrix[0,2]
        probat=pmatrix[0,3]
        probachild1=(probaa*node.children[1].sitearrays[siteindex][0]+
        probac*node.children[1].sitearrays[siteindex][1]+
        probag*node.children[1].sitearrays[siteindex][2]+
        probat*node.children[1].sitearrays[siteindex][3]) 
        #print 'probachild1=',probachild1
        
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild1=(probca*node.children[1].sitearrays[siteindex][0]+
        probcc*node.children[1].sitearrays[siteindex][1]+
        probcg*node.children[1].sitearrays[siteindex][2]+
        probct*node.children[1].sitearrays[siteindex][3]) 
        #print 'probcchild1=',probcchild1
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild1=(probga*node.children[1].sitearrays[siteindex][0]+
        probgc*node.children[1].sitearrays[siteindex][1]+
        probgg*node.children[1].sitearrays[siteindex][2]+
        probgt*node.children[1].sitearrays[siteindex][3])
        #print 'probgchild1=',probgchild1
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild1=(probta*node.children[1].sitearrays[siteindex][0]+
        probtc*node.children[1].sitearrays[siteindex][1]+
        probtg*node.children[1].sitearrays[siteindex][2]+
        probtt*node.children[1].sitearrays[siteindex][3])  
        #print 'probtchild1=',probtchild1

        """multiply probs for children 0 and 1 and get final array"""
       # print 'calculating array'
        array=(probachild0*probachild1,probcchild0*probcchild1,probgchild0*probgchild1,probtchild0*probtchild1)
        node.sitearrays=[]
        node.sitearrays.append(array)
        return array            

def findtips(node,tiplist=None):
    """takes a node as argument and returns a list containing all descendent nodes that do not yet have
    sitearrays associated to them. In more technical terms, i suppose that is to say that it finds the tips of
    the current tree after it has been pruned of the subtrees whose site prob arrays have alraedy been calculated"""
    """seems to be working perfectly!!!!"""
    """ANTICIPATED PROBLEM: WHEN I CALCULATE ARRAYS FOR THE SECOND, THIRD ETC SITES, NODES WILL ALREADY HAVE
    OBJECT SITEARRAY. ANTICIPATED SOLUTION: INSTEAD OF TESTING FOR EXISTENCE OF SITEARRAYS, TEST FOR EXISTENCE OF
    SITEARRAYS[SITEINDEX]"""
    if tiplist == None:
        tiplist=[]
    if len(node.children)==0 and hasattr(node,'sitearrays')==False: #if node is a (real) tip and does not yet have site array
        #print 'reached a tip:',node.name, 'Appending it to list'
        tiplist.append(node)
        #print 'list so far is:',tiplist
    elif hasattr(node,'sitearrays')==False and hasattr(node.children[0],'sitearrays')==True and hasattr(node.children[1],'sitearrays')==True: #second and third conditions assure that it is the last node without array
        #print 'reached a last node without sitearrays',node.name,'Appending it to list'
        tiplist.append(node)
        #print 'list so far is:',tiplist
    else:
        for child in node.children:
            findtips(child,tiplist=tiplist)
    return tiplist


def sitearraytree(node,siteindex):
    """takes a node (typically the root) root as argument and uses function sitearray (above) to calculate arrays for a given
    site at every node in tree. RECURSIVE"""
    """ the idea is to use my findtips function (above), that outputs a list of the last (closer to tips) nodes that do 
    not yet have site arrays, and then use my function sitearray to calculate arrays for nodes in that list. That is 
    implemented as a while loop that goes on until the root node itself has a sitearray (that is, hasattr(node,'sitearrays)
    bacomes True)."""
    """IT WORKS!!!!!"""
    while hasattr(node,'sitearrays')==False : 
        for tip in findtips(node) :
            sitearray(tip,siteindex)

      
        
    def treelike(root):
        """uses sitearraytree (above) to calculate arrays for every site in seq (i'll probably
        have to use smthg like a list of lists), then sums itens in each array
        to get site likelihoods, and finally multiplies these to get tree likelihood"""
        


