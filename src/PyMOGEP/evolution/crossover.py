# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''
import random

__all__ = ['crossoverPairs', 'crossoverOnePoint', 'crossoverTwoPoints', 'crossoverGene']

def crossoverPairs(popSize, crossoverRate):
    '''
    finding out which two chromosomes in the population should do crossover operation
    @param crossoverRate: crossover rate
    @return: ((p1a-index, p1b-index), (p2a-index, p2b-index), ...)
    '''
    if crossoverRate and popSize >= 2:
    # Choose and shuffle the individuals to participate in crossover
        orgs = [idx for idx in xrange(1, popSize) if random.random() < crossoverRate]
        random.shuffle(orgs)
            
        # Generate pairs of indexes.  If there is an odd man out, ignore him.
        for idx in xrange(0, len(orgs), 2):
            try:
                yield orgs[idx], orgs[idx+1]
            except IndexError:
                pass

def crossoverOnePoint(chro1, chro2):
        '''
        Produces two children via one-point crossover
        @param chro1: chro1 chromosome
        @param chro2: chro2 chromosome
        @return:      child 1, child 2
        '''
        genes1, genes2 = list(chro1.genes), list(chro2.genes)
        
        # Pick a gene and index to crossover at
        gene  = random.choice(xrange(len(genes1)))
        index = random.choice(xrange(len(genes1[gene])))
        
        # Construct new child genes
        child1 = genes1[gene].derive([(index, genes2[gene][index:])])
        child2 = genes2[gene].derive([(index, genes1[gene][index:])])
        genes1[gene], genes2[gene] = child1, child2
        return chro1._child(genes1), chro2._child(genes2)
    
def crossoverTwoPoints(chro1, chro2):
        '''
        Produces two children via two-point crossover
        @param chro1: chro1 chromosome
        @param chro2: chro2 chromosome
        @return:      child 1, child 2
        '''
        #total number of alleles in chro1 chromosome
        if len(chro1) < 2:
            return chro1, chro2

        genes1, genes2 = list(chro1.genes), list(chro2.genes)

        # Choose start and stop loci
        ind1, ind2 = random.sample(xrange(len(chro1)), 2)
        if ind1 > ind2:
            ind1, ind2 = ind2, ind1
        
        # Convert these to gene and allele numbers
        geneLength = len(chro1.genes[0])
        ind1, allele1 = divmod(ind1, geneLength)
        ind2, allele2 = divmod(ind2, geneLength)

        # Switch genes in between the modified genes
        if ind2 - ind1 > 1:
            start = ind1 + 1
            genes1[start:ind2], genes2[start:ind2] = \
                genes2[start:ind2], genes1[start:ind2]
        
        # And switch components of the start and stop genes
        child1 = genes1[ind1].derive([(allele1, genes2[ind1][allele1:])])
        child2 = genes2[ind1].derive([(allele1, genes1[ind1][allele1:])])
        genes1[ind1], genes2[ind1] = child1, child2
            
        child1 = genes1[ind2].derive([(0, genes2[ind2][:allele2])])
        child2 = genes2[ind2].derive([(0, genes1[ind2][:allele2])])
        genes1[ind2], genes2[ind2] = child1, child2
        
        return chro1._child(genes1), chro2._child(genes2)
    
def crossoverGene(chro1, chro2):
        '''
        Produces two children via full gene crossover
        @param chro1: chro1 chromosome
        @param chro2: chro2 chromosome
        @return:      child 1, child 2
        '''
        genes1, genes2 = list(chro1.genes), list(chro2.genes)

        # Choose a random gene
        gene = random.choice(xrange(len(genes1)))
        genes1[gene], genes2[gene] = genes2[gene], genes1[gene]
        return chro1._child(genes1), chro2._child(genes2)

