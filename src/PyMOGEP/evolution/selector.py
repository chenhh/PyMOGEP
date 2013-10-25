# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective chromosome is not able to use roulette wheel method,
using rank roulette wheeel or tournament
'''
import random, copy
from PyMOGEP.evolution.comparison import partialOrder

def uniformSelection(population):
    '''
    each chromosome has the same probability to be selected
    @param population, list of chromosome
    @return nextPopulation
    '''
    populationSize = len(population)
    return [copy.deepcopy(population[random.randint(0, populationSize-1)]) 
            for _ in xrange(populationSize)]
     

def binaryTournamentSelection(population):
    '''
    Two individuals are randomly chosen; 
    the fitter of the two is selected as a parent
    Note: utilizing fitnessesPartialOrder
    @param population, list of chromosome
    '''
#    print "tournament pop:", population
    populationSize = len(population)
    offSpring = [None] * populationSize 
    for idx in xrange(populationSize):
        jdx = random.randint(0, populationSize-1)
        kdx = random.randint(0, populationSize-1)
        #這邊必須使用copy, 否則nonDominatedSort會產生錯誤
        if partialOrder(population[jdx], population[kdx]):
            offSpring[idx] = copy.copy(population[jdx]) 
        else:
            offSpring[idx] = copy.copy(population[kdx])
              
    return offSpring