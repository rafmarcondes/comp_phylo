# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:21:17 2015

@author: Rafael
"""
from CTMCtree import sdd
import numpy
myq=numpy.matrix('-1.916 0.541 0.787 0.588; 0.148 -1.069 0.415 0.506; 0.286 0.170 -0.591 0.135; 0.525 0.236 0.594 -1.355') 


class ctmc(object):
    """         
    A class defining a continuous-time Markov chain (CTMC)
    """
          
    def __init__(self,q,staspa,v,simulations,times,nsites,starts,finals,sitemargprobs,seqmargprob):
        self.q=q 
        self.staspa=staspa
        self.v=v
        self.simulations=simulations
        self.times=times
        self.nsites=nsites
        self.starts=starts #starting states
        self.finals=finals #final states
        self.sitemargprobs=sitemargprobs
        self.seqmargprob=seqmargprob
        



    def simulate(self):
        for n in range(self.nsites) :
            """IMPORTANT:rows and columns in the q-matrix MUST be in the order ACGT
            IMPORTANT2: the q-matrix should be a numpy matrix object, a 'list of lists' won't work
            IMPORTANT3: this function is dependendt on my discrete sampling function, 'sdd', pasted above"""
            
            """ This method has no return statement. Rather, it simply updates 
            the ctmc attributes 'simulations', 'times', 'starts' and 'finals'"""
            import random
            statup=tuple(self.staspa) 
            """during the function, i'll have to remove elements from the 'states' list, but later reset the list 
            to its original composition. i'll use the 'statup' tuple as a 'backup'to do that"""
            chain=[]
            times=[]
            elapsedtime=0 #initialize a variable for the elapsed time
            if node.parent is None : # that is, if node is root
                currstate=random.choice(self.staspa) #CHANGE this TO SAMPLE INITIAL STATE FROM STATIONARY PROBS
            else : #IF NODE IS NOT ROOT 
                currstate=node.parent.seq[n]      
            chain.append(currstate)
            wt=0 #initialize the waiting time variable
            while elapsedtime<self.v :
                """this series of if statements chooses the appropriate lambda (rate of the 
                exponential distribution), given the currente state"""
                if currstate=='a' :
                    lambd=-(self.q.item(0,0))
                if currstate=='c' :
                    lambd=-(self.q.item(1,1))
                if currstate=='g' :
                    lambd=-(self.q.item(2,2))
                if currstate=='t' :
                    lambd=-(self.q.item(3,3))
                wt=random.expovariate(lambd) #use the appropriate lambda to draw a waiting time
                elapsedtime+=wt
                if elapsedtime>self.v : 
                    """if elapsedtime will exceed v with the addition of the last wt, 
                    then undo its addition to elapsed time and break the while loop, because v has been reached"""
                    elapsedtime-=wt
                    break
                times.append(wt)
                """this series of if/elif statements draws the next state"""       
                if currstate=='a' :
                    self.staspa.remove('a') #in order to not draw the same state as the current one
                    currstate=sdd(self.staspa,[self.q.item(1)/-self.q.item(0),self.q.item(2)/-self.q.item(0),self.q.item(3)/-self.q.item(0)])
                elif currstate=='c' :
                    self.staspa.remove('c')
                    currstate=sdd(self.staspa,[self.q.item(4)/-self.q.item(5),self.q.item(6)/-self.q.item(5),self.q.item(7)/-self.q.item(5)])
                elif currstate=='g' :
                    self.staspa.remove('g')
                    currstate=sdd(self.staspa,[self.q.item(8)/-self.q.item(10),self.q.item(9)/-self.q.item(10),self.q.item(11)/-self.q.item(11)])
                elif currstate=='t' :
                    self.staspa.remove('t')
                    currstate=sdd(self.staspa,[self.q.item(12)/-self.q.item(15),self.q.item(13)/-self.q.item(15),self.q.item(14)/-self.q.item(15)])
                chain.append(currstate)
                self.staspa=list(statup) #reset the 'states' list to its original composition
            times.append(self.v-elapsedtime) 
            """append to the times list the difference between the elapsed time so far and v, so that the 
            final elapsedtime, or sum(times), equals v"""
            self.starts+=(chain[0])
            self.finals+=(chain[-1])
            self.simulations.append(chain) #append the class attributes chains and times with the output of the simulations
            self.times.append(times)
    
    def margprobssites(self) :
        """outputs a list contining margprobs for each site in sequence"""
        sitemargprobs=[]
        import scipy
        pmatrix=scipy.linalg.expm(self.q*self.v)
        for i in range(self.nsites) :
            initial=self.starts[i]
            final=self.finals[i]
            iindex=self.staspa.index(initial)
            findex=self.staspa.index(final)
            mp=pmatrix[iindex,findex]
            sitemargprobs.append(mp)
        return sitemargprobs
        

        
    def margprobseq(self,v=None) :
        """margprob of whole seq, calculated by multiplying margprobs for each site.
                
        The keyword argument v is supposed to be used in situations where i don't want to use the self.v associated
        with the object and used to run the simulation, such as in calculating ML (method below), where the v used will
        be vCurr, vUp and vDown"""  
        if v is None :    
            sitemargprobs=[]
            for i in range(self.nsites) :
                initial=self.starts[i]
                final=self.finals[i]
                import scipy    
                pmatrix=scipy.linalg.expm(self.q*self.v)
                iindex=self.staspa.index(initial)
                findex=self.staspa.index(final)
                mp=pmatrix[iindex,findex]
                sitemargprobs.append(mp)
            seqmargprob=1
            for j in sitemargprobs:
                seqmargprob*=j
            return seqmargprob
        else:
            sitemargprobs=[]
            for i in range(self.nsites) :
                initial=self.starts[i]
                final=self.finals[i]
                import scipy    
                pmatrix=scipy.linalg.expm(self.q*v)
                iindex=self.staspa.index(initial)
                findex=self.staspa.index(final)
                mp=pmatrix[iindex,findex]
                sitemargprobs.append(mp)
            seqmargprob=1
            for j in sitemargprobs:
                seqmargprob*=j
            return seqmargprob
                

    def mlv(self,vStart,diff): 
        """vStart MUST be different from 0"""
        vCurr=vStart
        while diff>0.0001 :
            """CHANGE THE ABOVE TO BE AN ARGUMENT"""
            likeCurr=self.margprobseq(v=vCurr)
            vUp=vCurr+diff
            vDown=vCurr-diff
            if (vDown < 0):
                vDown = 0
            likeUp=self.margprobseq(v=vUp)
            likeDown=self.margprobseq(v=vDown)
            if likeDown>likeCurr :
                vCurr=vDown
                vUp=vCurr+diff
                vDown=vCurr-diff
            elif likeUp>likeCurr :
                vCurr=vUp
                vUp=vCurr+diff
                vDown=vCurr-diff
            else :
                diff=(diff/2.0)            
        return vCurr             

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
        node.markov=ctmc(staspa=['a','c','g','t'],q=myq,v=node.brl,simulations=[],times=[],nsites=30,starts='',finals='',sitemargprobs=[],seqmargprob=1)
        """ each node will have a ctmc object called 'markov' """        
        node.markov.simulate()
        node.seq=node.markov.finals
        """ in addition to 'node.markov.finals', seq in the node will also be stored in
        node.seq, for easier retrieval"""
        for child in node.children:
            self.simulate(child)
            
        """so far, the ctmc object attached to the node object seems to be working alright, as is running
        the simulation in the tree using recursion (the method above). What I still have to do, though, is
        somehow make the initial sequence in each node be the final sequence passed from the parent node. I think
        that will require changes in the ctmc.simulate method.
        """
       
        
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

print mytree.spA.seq

mytree.simulate(mytree.root)


print mytree.ancAB.seq
print mytree.spA.seq

print mytree.root.name
print mytree.spC.name

print mytree.root.children[0].name

print mytree.spB.parent.name
print mytree.root.parent

print mytree.spA.brl
