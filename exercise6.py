# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:21:17 2015

@author: Rafael
"""

from CTMCtree import ctmc 
import numpy
myq=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355') 

class Node(object):
   
    def __init__(self,markov=None,name="",parent=None,children=None,brl=0):
        self.name = name
        self.parent = None
        if children is None:
            self.children = []
        else:
            self.children = children
        self.brl=brl

        
            
    def printnames(self):
        """the method we physically acted in class, to print the names of all the tips
        descending from a node"""
        if len(self.children)==0:
            print self.name
        else:        
            for child in self.children:
                child.printnames()


class Tree:
    """
    Defines a class of phylogenetic tree, consisting of linked Node objects.
    """
    
    def __init__(self,root,spC,ancAB,spA,spB):
       
        self.root = Node(name="root") 
        self.spC = Node(name="SpeciesC",parent=self.root)
        self.root.children.append(self.spC)
        self.ancAB = Node(name="ancAB",parent=self.root)
        self.root.children.append(self.ancAB)
        self.spA = Node(name="SpeciesA",parent=self.ancAB)
        self.spB = Node(name="SpeciesB",parent=self.ancAB)
        self.ancAB.children.append(self.spA)
        self.ancAB.children.append(self.spB)
        # Now, let's add branch lengths to our Node objects (remember, these fields
        # can be added arbitrarily in Python). In the future, we should probably include
        # branch lengths in the Node constructor.
        self.spA.brl = 0.4
        self.spB.brl = 0.5
        self.spC.brl = 0.5
        self.ancAB.brl = 0.3
        self.root.brl = 0
        # We're also going to add lists to each node that will hold simulated
        # sequences.
        self.spA.seq = []
        self.spB.seq = []
        self.spC.seq = []
        self.ancAB.seq = []
        self.root.seq = []
        #self.setModels(self.root)

    
    def printNames(self,node):
        if len(node.children)==0:
            print node.name
        else:        
            for child in node.children:
                child.printNames()
 
   
    def treeLength(self,node):  # It works!!!
        """
        A method to calculate and return total tree length.
        """
        ttl=[]
        if not node.children: 
            return node.brl
        else:
            ttl.append(node.brl)
            for child in node.children:
                ttl.append(self.treeLength(child))
        return sum(ttl)


    #def newick(self,node):

    def setmodels(self,node,qmatrix):
        """
        This method of a Tree object defines a ctmc object associated with all
        nodes that have a branch length (i.e., all but the root).
        
        So far, the method also defines a ctmc associated with the root. This doesnt cause problems,
        because then the simulation of evolution in the root will just have the same initial and final
        sequences, given that its brl is 0. I guess this may be precious memory and time consumed, though, and
        plan on improving on this soon.
        """
        
        node.markov=ctmc(q=qmatrix,staspa=['a','c','g','t'],v=node.brl,simulations=[],times=[],nsites=30,starts='',finals='',sitemargprobs=[],seqmargprob=1)
        for child in node.children :
            self.setmodels(node=child,qmatrix=qmatrix)
        

    def treesimulate(self,node):
        """
        This method simulates evolution along the branches of a tree, taking
        the root node as its initial argument.
        """ 
        if node.parent==None :        
            node.markov.branchsimulate()
            node.seq=node.markov.finals
        else :
            node.markov.branchsimulate(parentnode=node.parent)
            node.seq=node.markov.finals
        """ in addition to 'node.markov.finals', seq in the node will also be stored in
        node.seq, for easier retrieval"""
        for child in node.children:
            self.treesimulate(child)
            
    def printseqs(self,node):
        """
        This method prints out the names of the tips and their associated
        sequences as an alignment (matrix).
        """
        if len(node.children)==0:
            print node.name, node.seq
        else:
            for child in node.children:
                self.printseqs(child)
                
                
            

        
"""
DEFINITIONS ABOVE, USES BELOW
"""


mytree=Tree('root','spC','ancAB','spA','spB')
mytree.spC.parent=mytree.root
mytree.ancAB.parent=mytree.root
mytree.spA.parent=mytree.ancAB
mytree.spB.parent=mytree.ancAB


mytree.root.printnames()
mytree.treeLength(mytree.root)

print mytree.root.name
print mytree.spC.name
print mytree.spA.brl
print mytree.spA.parent.name

for n in mytree.root.children:
    print n.name

print mytree.root.children[0].name

print mytree.spB.parent.name
print mytree.root.parent

print mytree.spA.brl
print mytree.root.brl


mytree.setmodels(node=mytree.root,qmatrix=myq)

mytree.treesimulate(mytree.root)

print mytree.root.seq
print mytree.ancAB.seq
print mytree.spA.seq



mytree.printseqs(mytree.root)
