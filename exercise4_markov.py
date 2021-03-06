"""
Exercise 4
Discrete-time Markov chains
@author: rafmarcondes
"""

"""
In this exercise, we will explore Markov chains that have discrete state spaces
and occur in discrete time steps. To set up a Markov chain, we first need to 
define the states that the chain can take over time, known as its state space.
To start, let's restrict ourselves to the case where our chain takes only two
states. We'll call them A and B.
"""

# Create a tuple that contains the names of the chain's states

staspa=("a","b")



"""
The behavior of the chain with respect to these states will be determined by 
the probabilities of taking state A or B, given that the chain is currently in 
A and B. Remember that these are called conditional probabilities (e.g., the 
probability of going to B, given that the chain is currently in state A is 
P(B|A).)

We record all of these probabilities in a transition matrix. Each row
of the matrix records the conditional probabilities of moving to the other
states, given that we're in the state associated with that row. In our example
row 1 will be A and row 2 will be B. So, row 1, column 1 is P(A|A); row 1, 
column 2 is P(B|A); row 2, column 1 is P(A|B); and row 2, column 2 is P(B|B). 
All of the probabilities in a ROW need to sum to 1 (i.e., the total probability
associated with all possibilities for the next step must sum to 1, conditional
on the chain's current state).

In Python, we often store matrices as "lists of lists". So, one list will be 
the container for the whole matrix and each element of that list will be 
another list corresponding to a row, like this: mat = [[r1c1,r1c2],[r2c1,r2c2]]. 
We can then access individual elements use two indices in a row. For instance,
mat[0][0] would return r1c1. Using just one index returns the whole row, like
this: mat[0] would return [r1c1,r1c2].

Define a transition matrix for your chain below. For now, keep the probabilties
moderate (between 0.2 and 0.8).
"""

# Define a transition probability matrix for the chain with states A and B

mymatrix=[[0.4,0.6],[0.5,0.5]]



# Try accessing a individual element or an individual row 
# Element

print(mymatrix[0][0]) #the element in row 0, column 0 (top left)

# Row

print(mymatrix[1]) #the top row


"""
Now, write a function that simulates the behavior of this chain over n time
steps. To do this, you'll need to return to our earlier exercise on drawing 
values from a discrete distribution. You'll need to be able to draw a random
number between 0 and 1 (built in to scipy), then use your discrete sampling 
function to draw one of your states based on this random number.
"""

# Import scipy U(0,1) random number generator

import random

# Paste or import your discrete sampling function


def sdd(events,probs):
    """This function is named after "Sample from Discrete Distribution". It takes as arguments two lists,
    'events' and 'probs', in which 'probs' contains floats 
    representing probabilities associated, respectively, with the itens in 'events', that can be strings, 
    floats, ints, or any data type. The two lists must be of same size."""
    
    import random
    nprobs=[x*1000 for x in probs] #so, here i multiply each float in 'probs' by 1000 and store the products in 'nprobs'
    newlist=[]
    for a in range(len(events)) : #then, in this loop, i create a list (newlist), in which each event appears 1000*its probability times
        b=nprobs[a]
        b=int(b)
        for c in range(b) :
            newlist.append(events[a]) 
    return (random.choice(newlist)) #and finally, i ramdonly sample an element from newlist


# Write your Markov chain simulator below. Record the states of your chain in 
# a list. Draw a random state to initiate the chain.



def markov(states,matrix,length) :
    
    """A function to simulate a discrete markov chain. Argument 'states' is a list or tuple representing the state space.
    'Matrix' is a transition probability matrix (a list of lists). 'length' is the number of steps you want in the chain.
    IMPORTANT: this function is dependent on my discrete sampling function, 'sdd'."""
    
    import random
    chain=[]
    currstate=random.choice(states)
    chain.append(currstate)
    for i in range ((length-1)) :
        for j in range (len(states)) :
            if currstate==states[j] :
                currstate=sdd(states,matrix[j])
                chain.append(currstate)
                break
    return(chain)
        

# Run a simulation of 10 steps and print the output.

print(markov(staspa,mymatrix,10))

# ----> Try to finish the above lines before Tues, Feb. 10th <----

# Now try running 100 simulations of 100 steps each. How often does the chain
# end in each state? How does this change as you change the transition matrix?


laststeps=[]
for n in range(100) : #run the 100 simulations
    chain=markov(("A","B"),[[0.1,0.9],[0.2,0.8]],100)  #assign each chain to 'chain'
    laststeps.append(chain[99]) #append the 10th (index 99) step of each chain to a laststeps list
print laststeps 

acount=0   
for i in laststeps : #count the number of times that the chains ended in "A"
    if i=="A" :
        acount+=1
print "OUT OF 100 CHAINS,", acount, "ENDED IN A"
    
"""It seems that the smaller I set the probability of transitioning to A ( A->A as well as B->A), the less
often the chains end in A"""


# Try defining a state space for nucleotides: A, C, G, and T. Now define a 
# transition matrix with equal probabilities of change between states.

nucs=("A","C","T","G")
nucmatrix=[[0.01,0.33,0.33,0.33],[0.1,0.3,0.3,0.3],[0.1,0.3,0.3,0.3],[0.1,0.3,0.3,0.3]]

# Again, run 100 simulations of 100 steps and look at the ending states. Then
# try changing the transition matrix.
         
laststeps=[]
for n in range(100) : #run the 100 simulations
    chain=markov(nucs,nucmatrix,100)  #assign each chain to 'chain'
    laststeps.append(chain[99]) #append the 10th (index 99) step of each chain to a laststeps list
print laststeps 

acount=0
tcount=0
ccount=0
gcount=0 
for i in laststeps : #count the number of times that the chains ended in "A"
    if i=="A" :
        acount+=1
    if i=="T" :
        tcount+=1 
    if i=="G" :
        gcount+=1
    if i=="C" :
        ccount+=1
print "OUT OF 100 CHAINS,", acount, "ENDED IN A"
print "OUT OF 100 CHAINS,", ccount, "ENDED IN C"
print "OUT OF 100 CHAINS,", tcount, "ENDED IN T"
print "OUT OF 100 CHAINS,", gcount, "ENDED IN G"

"""For an uniform transition probability matrix (each transition prob=0.25), each nucleotide tends to be
in the end of the chain 25% of times. As i reduce the probs of transitioning to a nucleotide, it appears less
often in the end of a chain"""
