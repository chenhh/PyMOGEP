# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective algorithm: NSGA-II

'''

from time import time
import numpy as np
from PyMOGEP.evolution.linker import defaultLinker
from PyMOGEP.sort import (DebSort, JensenSort)

class Population(object):
    '''population of GEP chromosomes'''
    # parameters of GEP
    
    mutationRate = 0.6 
    inversionRate = 0.1
    transpositionISRate = 0.1
    transpositionISLength = 1, 2, 3
    transpositionRISRate = 0.1
    transpositionRISLength = 1, 2, 3
    transpositionGeneRate = 0.1
    crossoverOnePointRate = 0.3
    crossoverTwoPointsRate = 0.3
    crossoverGeneRate = 0.1
    runStatistics = False
    
    df = None
    
    gen = property(lambda self: self._gen, doc='Generation number')
    bestFront = property(
        lambda self: self.ParetoFronts[0],
        doc='The best Pareto front Chromosome of the current generation'
    )

    def __init__(self, chromosome, popSize, headLength,
                 n_genes=1, linker=defaultLinker, verbose=False):
        '''
        @param chromosome, PyMOGEP.chromosome, user defined chromosome
        @param popSize, positive integer, population size
        @param headLength, positive integer, head length of a gene
        @param n_genes, positive integer, num. of genes of a chromsome
        @param linker, linker function for connecting genes
        Constructor
        '''
        assert popSize > 0 and headLength > 0 and n_genes > 0
        self.popSize = popSize
        self.headLength = headLength
        self.n_genes = n_genes
        self.linker = linker
        self.__age = 0
        self.selector = binaryTournamentSelection
        self.verbose = verbose
        
        # population initialization, 
        # note: we don't want the chromosomes with fitnesses are all zeros
        self.population = []
        t0 = time()
        while len(self.population) < popSize:
            chro = chromosome.randomChromosome(headLength, n_genes, linker)
            zero = True
            for val in chro.fitnesses:
                if val != 0.0:
                    zero = False
                    break
            if not zero:
                self.population.append(chro)
        
        if self.verbose:
            print "initialized random population, %.3f secs"%(time() - t0)
        
        # placeholder for next generation
        self._nextPopulation = [] 
        
        # placeholder for Pareto front
        self.n_objectives = self.population[0].n_objectives
        self.ParetoFronts = self._fastNonDominatedSort(self.population)
        [self._crowdingDistanceAssignment(front) for front in self.ParetoFronts]
        
        if self.verbose:
            print "Num. of objctives:", self.n_objectives
            print "initialize population, %.3f secs"%(time() - t0)

        # Compute stats about the initial generation
        if self.runStatistics:
            self.means = [float('inf')] * self.n_objectives
            self.stdevs = [float('inf')] * self.n_objectives
            self._updateStats()
      
    def _crowdingDistanceAssignment(self, nonDominatedSet):
        '''       
        The overall crowding-distance value is calculated as
        the sum of individual distance values corresponding to 
        each objective.
        @param nondominatedSet: one set of the Pareto front
        '''
        length = len(nonDominatedSet)
    
        # initialize distance
        for chro in nonDominatedSet:
            chro.crowdingDistance = 0.0
        
        # sorting by each objective(fitness)
        for idx in xrange(self.n_objectives):
            #sort idx-th objective (ascend)
            nonDominatedSet = sorted(nonDominatedSet,
                                key=lambda chro: chro.fitnesses[idx])
                           
            minObjValue = nonDominatedSet[0].fitnesses[idx]
            maxObjValue = nonDominatedSet[-1].fitnesses[idx]
            
            nonDominatedSet[0].crowdingDistance = float("inf")
            nonDominatedSet[-1].crowdingDistance = float("inf")
            
            for jdx in xrange(1, length - 1):
                if minObjValue != maxObjValue:
                #  in small pop, minObjVal may be the same as maxObjVal 
                    nonDominatedSet[jdx].crowdingDistance += (
                        nonDominatedSet[jdx + 1].fitnesses[idx] - 
                        nonDominatedSet[jdx - 1].fitnesses[idx]) / (
                        maxObjValue - minObjValue)
                else:
                    nonDominatedSet[jdx].crowdingDistance += (
                        nonDominatedSet[jdx + 1].fitnesses[idx] - 
                        nonDominatedSet[jdx - 1].fitnesses[idx])
                
        return nonDominatedSet
    
    def _fastNonDominatedSort(self, population):
        '''@return all Pareto fronts (list of non-dominated set)'''
        if self.n_objectives == 2:
            ParetoFronts = JensenSort.twoObjectivesSweepAlgorithm(population)
        else:
            ParetoFronts = JensenSort.highObjectivesNonDominatedSort(population)
        return ParetoFronts
    
    def evolution(self, population):
        '''
        -如果pop中有一個以上的chromosome的fitnesses為[0, 0] (2-obj)時，
        -將只保留一個chromosome, 再用Pareto front rank 1中的chromosome
        -來演化出fitnesses非[0,0]的chromosome
        '''
        
        population = [invert(chro, self.inversionRate) for chro in population]
        
        # Insertion Sequence transposition
        population = [transposeIS(chro, random.choice(self.transpositionISLength),
                        self.transpositionISRate) for chro in population ]
        
        # Root Insert Sequence transposition
        population = [transposeRIS(chro, random.choice(self.transpositionRISLength),
                        self.transpositionRISRate) for chro in population]
        # Gene transposition
        population = [transposeGene(chro, self.transpositionGeneRate) for chro in population]
           
        # mutation
        population = [mutation(chro, self.mutationRate) 
                      for chro in population if self.mutationRate]
        
        # crossover
        for idx, jdx in crossoverPairs(self.popSize, self.crossoverOnePointRate):
            par1, par2 = population[idx], population[jdx]            
            child1, child2 = crossoverOnePoint(par1, par2)
            population[idx], population[jdx] = child1, child2

        for idx, jdx in crossoverPairs(self.popSize, self.crossoverTwoPointsRate):
            par1, par2 = population[idx], population[jdx]            
            child1, child2 = crossoverTwoPoints(par1, par2)
            population[idx], population[jdx] = child1, child2

        for idx, jdx in crossoverPairs(self.popSize, self.crossoverGeneRate):
            par1, par2 = population[idx], population[jdx]            
            child1, child2 = crossoverGene(par1, par2)
            population[idx], population[jdx] = child1, child2
        
        return population
    
    def evolveNSGAII(self):
        # produce offspring        
        currentOffspring = self.selector(self.population)
        
        currentOffspring = self.evolution(currentOffspring)    
        
        mixedPopulation = self.population + currentOffspring
        mixedParetoFronts = self._fastNonDominatedSort(mixedPopulation)
        
        assert sum(len(front) for front in mixedParetoFronts) == 2 * self.popSize
        
        # fill out the next population set
        self._nextPopulation = []
        idx = 0 
        while (len(self._nextPopulation) + len(mixedParetoFronts[idx])) <= self.popSize:
            self._nextPopulation.extend(self._crowdingDistanceAssignment(mixedParetoFronts[idx]))
            idx += 1
            
        # 補齊chromosome to next popluation
        self._crowdingDistanceAssignment(mixedParetoFronts[idx])

        # descending order
        mixedParetoFronts[idx].sort(cmp=partialOrder, reverse=True)
        reminderLength = self.popSize - len(self._nextPopulation)
        self._nextPopulation.extend(mixedParetoFronts[idx][:reminderLength])
        
        assert len(self._nextPopulation) == self.popSize

        # update population, Pareto fronts and crowding distance
        self._nextPopulation, self.population = [], self._nextPopulation
        self.ParetoFronts = self._fastNonDominatedSort(self.population)
        [self._crowdingDistanceAssignment(front) for front in self.ParetoFronts]
        
        # update information
        self.__age += 1
        if self.runStatistics:
            self._updateStats()
 
    def solve(self, generations):
        '''
        Cycles a number of generations. Stops if self.solved()
        @param generations: # of genrations to give up after
        '''
        for _ in xrange(generations):
            t1 = time()
            self.evolveNSGAII()
            if len(self.bestFront) > 0 and all([front.solved for front in self.bestFront]):
                break
            print self
            print "Generation:[{0}] Best Pareto front:".format(self.gen)
            for chro in self.bestFront:
                print  chro.fitnesses, chro
            print "best front size:{0}, execution time:{1:.3f} secs".format(
                len(self.bestFront), time() - t1)
        
