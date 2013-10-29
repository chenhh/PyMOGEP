# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective algorithm: NSGA-II

M. T. Jensen, "Reducing the run-time complexity of multiobjective EAs: 
The NSGA-II and other algorithms," Evolutionary Computation, IEEE Transactions on, 
vol. 7, pp. 503-515, 2003.
'''

import bisect

def twoObjectivCmpFunc(chro1, chro2):
    '''
    if chromosome1.fitnesses dominating chromosome2.fitnesses, then return 1
    elif chromosome1.fitnesses == chromosome2.fitnesses, then return 0
    else return -1
    '''
    assert chro1.n_objectives == chro2.n_objectives == 2
    
    if (chro1.fitnesses[0] < chro2.fitnesses[0]) or (
     chro1.fitnesses[0] == chro2.fitnesses[0] and 
     chro1.fitnesses[1] < chro2.fitnesses[1]):
        return 1
    elif chro1.fitnesses == chro2.fitnesses:
        return 0
    else:
        return -1


def twoObjectivesSweepAlgorithm(population):
    '''
    @param population, PyMOGEP.population
    time complexity O(NlogN), N: population size
    '''
    
    #decreasing (dominating), O(NlogN)
    population.sort(cmp=twoObjectivCmpFunc, reverse=True)
    ParetoFronts =[[population[0]],]
    
    frontCounter = 0
    for chro in population[1:]:
        #check if current chromosome is dominated by 
        #any chromosome in ParetoFronts[frontCounter]
        isDominated = False
        for frontChro in ParetoFronts[frontCounter]:
            if frontChro.dominating(chro):
                isDominated = True
                break
            
        if isDominated:
            #because the chromosome is isDominated by one chromosome
            #in the Pareto front level i, the chromosome must at least 
            #belong to Pareto front level (i+1)
            ParetoFronts.append([])
            frontCounter += 1
            ParetoFronts[frontCounter].append(chro)
            
        else:
            #find lowest front b such that
            # any chromosome in Paretofronts[b] not dominating the chromosome
            b = 0
            for level in xrange(frontCounter+1):
                if not any(chro2.dominating(chro) for chro2 in ParetoFronts[level]):
                    b = level
                    break
            ParetoFronts[b].append(chro)
    return ParetoFronts


def highObjectivesNonDominatedSort(population):
    '''@param population, PyMOGEP.population'''
    
    assert population[0].n_objectives > 2
    
    #initialize rank
    for chro in population:
        chro.ParetoRank = 1
    
    ND_helper_A(population)
    
    #allocating each point to corresponding front    
    n_front = max(chro.ParetoRank for chro in population)
    ParetoFronts = [ [] for _ in xrange(n_front)]
    [ParetoFronts[chro.ParetoRank-1].append(chro) for chro in population]
    
    return ParetoFronts


def splitSet(population):
    '''
    split population according to m-th objective value.
    @return median index of m-th objective value.
    '''
    popSize = len(population)
    n_objectives = population[0].n_objectives
    
    midPopIdx  = (popSize-1)/2 if popSize % 2 else popSize/2 - 1

    #all object in set H will not dominating L
    for objectiveIdx in reversed(xrange(n_objectives)):
        population.sort(key = lambda chro: chro.fitnesses[objectiveIdx])
        values = [chro.fitnesses[objectiveIdx] for chro in population]
        
        medianIdx = bisect.bisect_right(values, values[midPopIdx])
        if medianIdx != popSize:
            return medianIdx-1, objectiveIdx
        else:
            #we can't split the list after index midPopIdx, 
            #hence, we try to split it beroe the midPopIdx
            medianIdx = bisect.bisect_left(values, values[midPopIdx])
            if medianIdx != 0:
                return medianIdx-1, objectiveIdx
    
    #all fitness values in the population are the same
    return -1, -1


def ND_helper_A(population):
    '''
    split population to two sets, L and H, according to objectives[objectiveIdx]
    1. find the median of objectives[objectiveIdx]
    2. the objectives[objectiveIdx] of all chromosomes in set L are less than
       the median 
    3. the objectives[objectiveIdx] of all chromosomes in set H are large than
       the median
       
    the chromosome in set H are impossible to dominate the chromosome in set L
    (because the m-th objective value of chromosome in set H are larger than 
    the m-th objectivce value of chromsome in set L)
    '''
    popSize = len(population)
    
    if popSize == 2:
        # stop condition
        if population[0].dominating(population[1]):
            population[1].ParetoRank = max(population[0].ParetoRank+1, 
                                           population[1].ParetoRank)
        elif  population[1].dominating(population[0]):
            population[0].ParetoRank = max(population[0].ParetoRank, 
                                           population[1].ParetoRank+1)
    
    elif popSize > 2:
        medianIdx, objectiveIdx = splitSet(population)
        if medianIdx == -1:
            #all fitess values in the population are the same
            return
        
        L, H = population[:medianIdx+1], population[medianIdx+1:]
        
        #check dominating relations in set L, and computing Pareto rank
        ND_helper_A(L)
        #merge
        ND_helper_B(L, H, objectiveIdx-1)
        #check dominating relations in set H
        ND_helper_A(H)
        

def ND_helper_B(L, H, objectiveIdx):
    '''
    the procedure assigns Pareto rank to the solutions in H 
    according to the solutions in L.
    assumpting all Pareto rank of the chromosome in set L are assigned.
    
    '''
    if len(L) == 1:
        #stop condition
        for chro in H:
            if L[0].dominating(chro):
                chro.ParetoRank = max(chro.ParetoRank, L[0].ParetoRank+1)
                
    elif len(H) == 1:
        #stop condition 
        for chro in L:
            if chro.dominating(H[0]):
                H[0].ParetoRank = max(chro.ParetoRank+1, H[0].ParetoRank)
                
    elif 0 <= objectiveIdx < L[0].n_objectives:
        #for better performance, we don't split the set deeply. 
        L.sort(key=lambda chro: chro.ParetoRank, reverse=True)  #descending
        for chroH in H:
            for chroL in L:
                if chroL.dominating(chroH):
                    chroH.ParetoRank = max(chroL.ParetoRank+1, chroH.ParetoRank)
                    break
  