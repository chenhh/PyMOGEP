# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

return type of each method is list of the following format:
(alleleIdx, [new allele1,new allele2,...])
'''
import random

__all__ = ['crossoverPairs', 'crossoverOnePoint', 
           'crossoverTwoPoints', 'crossoverGene']

def crossoverPairs(popSize, crossoverRate):
    '''
    finding out which two chromosomes in the population should do 
    crossover operation
    
    @param crossoverRate: positive float
    @return: ((pair1-a-idx, pair1-b-idx), 
              (pair2-a-idx, pair2-b-idx), ...)
    '''
    assert 0 < crossoverRate <= 1.
    
    if crossoverRate and popSize >= 2:
        indices = [idx for idx in xrange(popSize) 
                   if random.random() < crossoverRate]
        random.shuffle(indices)
        
        if len(indices) % 2: 
            indices = indices[:-1]
        
        for idx in xrange(0, len(indices), 2):
            yield indices[idx], indices[idx+1]
          

def crossoverOnePoint(chro1, chro2):
        '''
        Produces two children via one-point crossover
        it will not produce illegel genes.
        
        @param chro1, PyMOGEP.chromosome
        @param chro2:  PyMOGEP.chromosome
        @return: child 1, child 2
        '''
        genes1, genes2 = list(chro1.genes), list(chro2.genes)
        
        # Pick a geneIdx and alleleIdx for crossover
        geneIdx  = random.choice(xrange(len(genes1)))
        alleleIdx = random.choice(xrange(len(genes1[geneIdx])))
        
        # Construct new child genes
        child1 = genes1[geneIdx].modify([(alleleIdx, genes2[geneIdx][alleleIdx:])])
        child2 = genes2[geneIdx].modify([(alleleIdx, genes1[geneIdx][alleleIdx:])])
        genes1[geneIdx], genes2[geneIdx] = child1, child2
        
        return chro1.newInstance(genes1), chro2.newInstance(genes2)


def crossoverTwoPoints(chro1, chro2):
        '''
        Produces two children via two-point crossover
        it will not produce illegel genes.
        
        @param chro1: PyMOGEP.chromosome
        @param chro2: PyMOGEP.chromosome
        @return: child 1, child 2
        '''
        #total number of alleles in chro1 chromosome
        if len(chro1) < 2:
            return chro1, chro2

        genes1, genes2 = list(chro1.genes), list(chro2.genes)

        #start and stop index（on chromosome)
        idx1, idx2 = random.sample(xrange(len(chro1)), 2)
        if idx1 > idx2:
            idx1, idx2 = idx2, idx1
        
        # Convert these to gene and allele numbers
        geneIdx1, alleleIdx1 = divmod(idx1, len(chro1.genes[0]))
        geneIdx2, alleleIdx2 = divmod(idx2, len(chro1.genes[0]))

        # Switch genes in between the modified genes
        if idx2 - idx1 > 1:
            start = idx1 + 1
            genes1[start:idx2], genes2[start:idx2] = \
                genes2[start:idx2], genes1[start:idx2]
        
        #switch alleles of the start and stop genes
        child1 = genes1[geneIdx1].modify([(alleleIdx1, genes2[geneIdx1][alleleIdx1:])])
        child2 = genes2[geneIdx1].modify([(alleleIdx1, genes1[geneIdx1][alleleIdx1:])])
        genes1[geneIdx1], genes2[geneIdx1] = child1, child2
            
        child1 = genes1[geneIdx2].modify([(0, genes2[geneIdx2][:alleleIdx2])])
        child2 = genes2[geneIdx2].modify([(0, genes1[geneIdx2][:alleleIdx2])])
        genes1[geneIdx2], genes2[geneIdx2] = child1, child2
        
        return chro1.newInstance(genes1), chro2.newInstance(genes2)


def crossoverGene(chro1, chro2):
        '''
        Produces two children via full geneIdx crossover
        it will not produce illegel genes.
        
        @param chro1: PyMOGEP.chromosome
        @param chro2: PyMOGEP.chromosome
        @return: child 1, child 2
        '''
        genes1, genes2 = list(chro1.genes), list(chro2.genes)

        # Choose a random geneIdx
        geneIdx = random.choice(xrange(len(genes1)))
        genes1[geneIdx], genes2[geneIdx] = genes2[geneIdx], genes1[geneIdx]
        return chro1.newInstance(genes1), chro2.newInstance(genes2)

