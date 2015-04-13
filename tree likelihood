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
    """takes as arguments a node and a site index (of seq in that node), and outputs a prob array 
    (in order a,c,g,t) for that site, based on child nodes"""
    """an attribute (a list) of node object is created, called sitearrays, that will contain prob
    arrays for each site in that node (thus it will be a list of lists)"""  
    """this function seems to be working perfectly"""
    if len(node.children)==0: #if node is a tip
        print 'len(children) is 0, node is a tip'
        if node.seq[siteindex]=='a':
            print 'nc at site at tip is a'
            array=(1,0,0,0)
        if node.seq[siteindex]=='c':
            print 'nc at site at tip is c'
            array=(0,1,0,0)
        if node.seq[siteindex]=='g':
            print 'nc at site at tip is g'
            array=(0,0,1,0)
        if node.seq[siteindex]=='t':
            print 'nc at site at tip is t'
            array=(0,0,0,1)
        node.sitearrays=[]
        node.sitearrays.append(array)
        return array
    else: #if node is not a tip
        """node is not a tip. calculations for child 0: """
        print 'starting calculatiosn for child 0'
        pmatrix=scipy.linalg.expm(node.children[0].markov.q*node.children[0].brl)#get pmatrix based on child's brl
        print 'pmatrix for child 0', pmatrix
        
        probaa=pmatrix[0,0]#prob of a->a transition along branch linking current node to child 0
        probac=pmatrix[0,1]#prob of a->c transition
        probag=pmatrix[0,2]# etc
        probat=pmatrix[0,3]
        probachild0=(probaa*node.children[0].sitearrays[siteindex][0]+ #multiplies tranition probs by prob that child has respective nc at site, and sums up to get prob of a at this site at this node, considering only child 0
        probac*node.children[0].sitearrays[siteindex][1]+
        probag*node.children[0].sitearrays[siteindex][2]+
        probat*node.children[0].sitearrays[siteindex][3]) 
        print 'probachild0=', probachild0
    
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild0=(probca*node.children[0].sitearrays[siteindex][0]+
        probcc*node.children[0].sitearrays[siteindex][1]+
        probcg*node.children[0].sitearrays[siteindex][2]+
        probct*node.children[0].sitearrays[siteindex][3]) 
        print 'probcchild0=',probcchild0
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild0=(probga*node.children[0].sitearrays[siteindex][0]+
        probgc*node.children[0].sitearrays[siteindex][1]+
        probgg*node.children[0].sitearrays[siteindex][2]+
        probgt*node.children[0].sitearrays[siteindex][3])
        print 'probgchild0=',probgchild0
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild0=(probta*node.children[0].sitearrays[siteindex][0]+
        probtc*node.children[0].sitearrays[siteindex][1]+
        probtg*node.children[0].sitearrays[siteindex][2]+
        probtt*node.children[0].sitearrays[siteindex][3])  
        print 'probtchild0=',probtchild0
                      
        """calculations for child 1:"""
        print 'starting calculatiosn for child 1'
        pmatrix=scipy.linalg.expm(node.children[1].markov.q*node.children[1].brl)#get pmatrix based on child's brl
        print 'pmatrix for child 0', pmatrix
        
        probaa=pmatrix[0,0]
        probac=pmatrix[0,1]
        probag=pmatrix[0,2]
        probat=pmatrix[0,3]
        probachild1=(probaa*node.children[1].sitearrays[siteindex][0]+
        probac*node.children[1].sitearrays[siteindex][1]+
        probag*node.children[1].sitearrays[siteindex][2]+
        probat*node.children[1].sitearrays[siteindex][3]) 
        print 'probachild1=',probachild1
        
        probca=pmatrix[1,0]
        probcc=pmatrix[1,1]
        probcg=pmatrix[1,2]
        probct=pmatrix[1,3]
        probcchild1=(probca*node.children[1].sitearrays[siteindex][0]+
        probcc*node.children[1].sitearrays[siteindex][1]+
        probcg*node.children[1].sitearrays[siteindex][2]+
        probct*node.children[1].sitearrays[siteindex][3]) 
        print 'probcchild1=',probcchild1
    
        probga=pmatrix[2,0]
        probgc=pmatrix[2,1]
        probgg=pmatrix[2,2]
        probgt=pmatrix[2,3]
        probgchild1=(probga*node.children[1].sitearrays[siteindex][0]+
        probgc*node.children[1].sitearrays[siteindex][1]+
        probgg*node.children[1].sitearrays[siteindex][2]+
        probgt*node.children[1].sitearrays[siteindex][3])
        print 'probgchild1=',probgchild1
    
        probta=pmatrix[3,0]
        probtc=pmatrix[3,1]
        probtg=pmatrix[3,2]
        probtt=pmatrix[3,3]
        probtchild1=(probta*node.children[1].sitearrays[siteindex][0]+
        probtc*node.children[1].sitearrays[siteindex][1]+
        probtg*node.children[1].sitearrays[siteindex][2]+
        probtt*node.children[1].sitearrays[siteindex][3])  
        print 'probtchild1=',probtchild1

        """multiply probs for children 0 and 1 and get final array"""
        print 'calculating array'
        array=(probachild0*probachild1,probcchild0*probcchild1,probgchild0*probgchild1,probtchild0*probtchild1)
        node.sitearrays=[]
        node.sitearrays.append(array)
        return array            

def findtips(node,tiplist=None):
    """takes a node as argument and returns a list containing all tips descending from that node"""
    if tiplist == None:
        tiplist=[]
    if len(node.children)==0: #if node is tip
        #print 'reached a tip:',node.name
        tiplist.append(node)
        #print 'tip list so far is',tiplist
    else:
        #print 'just tried',node.name,'and it is not tip'
        for child in node.children:
            findtips(child,tiplist=tiplist)
    return tiplist


def sitearraytree(node,siteindex,recursing=False):
    """takes root as argument and uses function sitearray (above) to calculate arrays for a given
    site at every node in tree. RECURSIVE"""
    """first, i have to go up the tree until i find a tip (a node without children).
    when i find this, use sitearray."""
    """i'm trying to do that by first creating a list of all tips in the tree,
    by using my findtips function"""
    """is this gonna work? may be better to somehow find all the nodes descending from node that not yet have arrays,
    instead of finding tips"""
    """yeah... not working... try that approach"""
    tips=findtips(node)
    for tip in tips:
        if len(tip.sitearrays)==0: #this is not working!!
            tips.remove(tip) #if site array has already been calculated for that node, remove it from this list (whose members will have arrays calculated below)
    for tip in tips:
        sitearray(tip,siteindex)
        sitearraytree(tip,siteindex)
    

    

        
    def treelike(root):
        """uses sitearraytree (above) to calculate arrays for every site in seq (i'll probably
        have to use smthg like a list of lists), then sums itens in each array
        to get site likelihoods, and finally multiplies these to get tree likelihood"""
        




#TESTS BELOW
    
myq=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355')


mytree=Tree('root','spC','ancAB','spA','spB')
mytree.spC.parent=mytree.root
mytree.ancAB.parent=mytree.root
mytree.spA.parent=mytree.ancAB
mytree.spB.parent=mytree.ancAB


mytree.spA.seq='aaacgtgcacgtgcattgtctgacggtcac'
mytree.spB.seq='aaacgtgcacgtgcattgcgacggggtcac'
mytree.spC.seq='aaacgcagtcgtgcattttgacggggtcac'




mytree.setmodels(mytree.root,myq)

#mytree.treesimulate(mytree.root)

print mytree.ancAB.seq
print mytree.ancAB.markov.starts
mytree.ancAB.markov.margprobssites()
print mytree.ancAB.markov.sitemargprobs
mytree.ancAB.markov.margprobseq()

mytree.printseqs(mytree.root)

print len(mytree.root.children)
    
#mylike=likelihood(tree=mytree,tips=[])

#print mylike



sitearray(mytree.spC,0)
sitearray(mytree.spA,0)
sitearray(mytree.spB,0)
print mytree.spC.array

sitearray(mytree.ancAB,0)  


sitearray(mytree.root,0)  

print mytree.ancAB.sitearrays[0]

sitearraytree(mytree.root,0)
    
