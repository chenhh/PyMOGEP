# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective algorithm: NSGA-II

'''

from time import time
import random
from PyMOGEP.evolution.linker import defaultLinker
from PyMOGEP.sort import (JensenSort, DebSort)
from PyMOGEP.evolution.selector import binaryTournamentSelection
from PyMOGEP.evolution.crossover import *
from PyMOGEP.evolution.mutator import *
from PyMOGEP.evolution.transposer import *
from PyMOGEP.evolution.comparison import *


class Population(object):
    '''population of GEP chromosomes'''
    
    # parameters of GEP
    mutationRate = 0.6 
    inversionRate = 0.1
    transISRate = 0.1
    transISLength = 1, 2, 3
    transRISRate = 0.1
    transRISLength = 1, 2, 3
    transGeneRate = 0.1
    crossoverOnePointRate = 0.3
    crossoverTwoPointsRate = 0.3
    crossoverGeneRate = 0.1
    df = None   #data frame
    
    gen = property(lambda self: self._gen, doc='Generation')
    bestFront = property(
        lambda self: self.ParetoFronts[0],
        doc='The best Pareto front'
    )

    def __init__(self, chro, popSize, headLength, n_genes=1, n_elites=1,
                 linker=defaultLinker, RNCGenerator=None, verbose=False):
        '''
        @param chro, PyMOGEP.chromosome, user defined chromosome
        @param popSize, positive integer, population size
        @param headLength, positive integer, head length of a gene
        @param n_genes, positive integer, number of gene of a chromsome
        @param n_elites, positive integer, number of better front preserved 
                                           in each generation.
        @param linker, PyMOGEP.evolution.linker, 
                       linker function for connecting genes
        @param RNCGenerator, random number generator for RNC algorithm
        '''
        assert popSize > 0 and headLength > 0 and n_genes > 0
        self.popSize = popSize
        self.headLength = headLength
        self.n_genes = n_genes
        self.n_elites = n_elites
        self.linker = linker
        self._gen = 0
        self.selector = binaryTournamentSelection
        self.RNCGenerator = RNCGenerator
        self.verbose = verbose
        
        #population initialization, each chromosome with different fitness values.
        self.population = []
        t0 = time()
        fitness_set = set()
        while len(self.population) < popSize:
            chro = chro.randomChromosome(headLength, n_genes, linker, self.RNCGenerator)
            if chro.fitnesses not in fitness_set:
                fitness_set.add(chro.fitnesses)
                self.population.append(chro)
        
        if self.verbose:
            for chro in self.population:
                print "ID:[%s], fitnesses:%s"%(chro.chromosomeID, chro.fitnesses)
                print chro
            print "random population initialized, %.3f secs"%(time() - t0)
        
        # placeholder for next generation
        self._nextPopulation = [] 
        
        # placeholder for Pareto front
        self.n_objectives = self.population[0].n_objectives
        self.ParetoFronts = self._fastNonDominatedSort(self.population)
        [self._crowdingDistanceAssignment(front) for front in self.ParetoFronts]
        
        if self.verbose:
            print "n_objctives:", self.n_objectives
            print "initialize population, %.3f secs"%(time() - t0)


    def _crowdingDistanceAssignment(self, nonDominatedSet):
        '''       
        The overall crowding-distance value is calculated as
        the sum of individual distance values corresponding to 
        each objective.
        the larger distance, the better chromosome
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
        if self.n_objectives == 1:
            population.sort(key=lambda chro: chro.fitnesses[0]) #ascending
            ParetoFronts = [[population[0]], ]
            population[0].ParetoRank = 1
            for idx, chro in enumerate(population[1:]):
                if chro.fitnesses[0] == population[idx].fitnesses[0]:
                    chro.ParetoRank = population[idx].ParetoRank
                    ParetoFronts[chro.ParetoRank-1].append(chro)
                else:
                    ParetoFronts.append([chro,])
                    chro.ParetoRank = population[idx].ParetoRank + 1
               
        if self.n_objectives == 2:
            ParetoFronts = JensenSort.twoObjectivesSweepAlgorithm(population)
        elif self.n_objectives > 2:
            ParetoFronts = JensenSort.highObjectivesNonDominatedSort(population)

        return ParetoFronts
    
    
    def evolution(self, population):
        #inversion
        population = [inversion(chro, self.inversionRate) for chro in population]
        
        # Insertion Sequence transposition
        population = [transposeIS(chro, random.choice(self.transISLength),
                        self.transISRate) for chro in population ]
        
        # Root Insert Sequence transposition
        population = [transposeRIS(chro, random.choice(self.transRISLength),
                        self.transRISRate) for chro in population]
        
        # Gene transposition
        population = [transposeGene(chro, self.transGeneRate) for chro in population]
           
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
    
    
    def evolve(self):
        '''execute the following procedure in each generation'''
        # produce offspring        
        offspring = self.selector(self.population)
        offspring = self.evolution(offspring)    
        
        mixedPopulation = self.population + offspring
        mixedParetoFronts = self._fastNonDominatedSort(mixedPopulation)
        
        assert sum(len(front) for front in mixedParetoFronts) == 2 * self.popSize
        
        # fill out the next population set
        self._nextPopulation = []
        
        #preserve elite Pareto front
        for idx in xrange(self.n_elites):
            self._nextPopulation.extend(
                self._crowdingDistanceAssignment(mixedParetoFronts[idx])
            )
        
        idx = self.n_elites 
        while (len(self._nextPopulation) + len(mixedParetoFronts[idx])) <= self.popSize:
            self._nextPopulation.extend(
                    self._crowdingDistanceAssignment(mixedParetoFronts[idx])
                    )
            idx += 1
        self._crowdingDistanceAssignment(mixedParetoFronts[idx])

        if len(self._nextPopulation) > self.popSize:
            self._nextPopulation = self._nextPopulation[:self.popSize]
        
        if len(self._nextPopulation) < self.popSize:
            #fullfill the popSize
            mixedParetoFronts[idx].sort(cmp=partialOrder, reverse=True)
            reminderLength = self.popSize - len(self._nextPopulation)
            self._nextPopulation.extend(mixedParetoFronts[idx][:reminderLength])
        
        assert len(self._nextPopulation) == self.popSize

        # update population, Pareto fronts and crowding distance
        self._nextPopulation, self.population = [], self._nextPopulation
        self.ParetoFronts = self._fastNonDominatedSort(self.population)
        [self._crowdingDistanceAssignment(front) for front in self.ParetoFronts]
        
        # update information
        self._gen += 1
     
 
    def solve(self, n_generation):
        '''
        execute a number of generations. Stops if self.solved()
        @param generation: number of generation for evolving
        '''
        for _ in xrange(n_generation):
            t0 = time()
            self.evolve()
            
            if len(self.bestFront) > 0 and all(chro.solved for chro in self.bestFront):
                break
            
            print "Generation[%s] %.3f secs, Best Pareto front size:%s"%(
                    self.gen, time()-t0, len(self.bestFront))
            
            if self.verbose:
                for chro in self.bestFront:
                    print  "1st rank:", chro.fitnesses
        
