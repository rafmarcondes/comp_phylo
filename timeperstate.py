def timeperstate(chainandtimes) :
    """a function that calculates the times spent in each state, given a markov chain and the waiting times
    between transitions.
    INPUT: a single argument, which must be a list or tuple in the format (chain,times), in which chain and times
    also are lists or tuples. This is the format of the output of my contmark function
    OUTPUT: a list with 4 itens, each being the total time spent in states a, c, g and t, in that order"""
    chain=chainandtimes[0]
    times=chainandtimes[1]
    ttime=0
    atime=0
    ctime=0
    gtime=0
    for n in range(0,len(chain)):
            if chain[n]=='a' :
                atime+=times[n]
            if chain[n]=='c' :
                ctime+=times[n]
            if chain[n]=='t' :
                ttime+=times[n]
            if chain[n]=='g' :
                gtime+=times[n]
    actgtimes=[]
    actgtimes.append(atime)
    actgtimes.append(ctime)
    actgtimes.append(gtime)
    actgtimes.append(ttime)
    return actgtimes
