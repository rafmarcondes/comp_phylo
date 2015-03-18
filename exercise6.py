# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:21:17 2015

@author: Rafael
"""
from Exercise5_CTMC import sdd
from Exercise5_CTMC import ctmc
import numpy
myq=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355') 


class Node(object):
   
    def __init__(self,name="",parent=None,children=None,brl=0):
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
       
        self.root = Node("root") 
        self.spC = Node("SpeciesC",parent=self.root)
        self.root.children.append(self.spC)
        self.ancAB = Node("ancAB",parent=self.root)
        self.root.children.append(self.ancAB)
        self.spA = Node("SpeciesA",parent=self.ancAB)
        self.spB = Node("SpeciesB",parent=self.ancAB)
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
                child.printnames()
 
   
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


    def simulate(self,node):
        """
        This method simulates evolution along the branches of a tree, taking
        the root node as its initial argument.
        """
        node.evolution=ctmc(staspa=['a','c','g','t'],q=myq,v=node.brl,simulations=[],times=[],nsites=30,starts='',finals='',sitemargprobs=[],seqmargprob=1)
        """ each node will have a ctmc object called 'evolution' """        
        node.evolution.simulate()
        node.seq=node.evolution.finals
        """ in addition to 'node.evolution.finals', the sequence in the node will also be stored in
        node.seq, for easier retrieval"""
        for child in node.children:
            self.simulate(child)
            
        """so far, the ctmc object attached to the node objects seems to be working alright, as is running
        the simulation in the tree using recursion (the method above). What I still have to do, though, is
        somehow make the initial sequence in each node be the final sequence passed from the parent node. I think
        that will require changes in the ctmc.simulate method.
        """
       
        
"""
DEFINITIONS ABOVE, USES BELOW
"""




mytree=Tree('root','spC','ancAB','spA','spB')

mytree.root.printnames()
mytree.treeLength(mytree.root)

print mytree.spA.seq

mytree.simulate(mytree.ancAB)

print mytree.ancAB.evolution.finals
print mytree.ancAB.seq



print mytree.spA.evolution.finals
