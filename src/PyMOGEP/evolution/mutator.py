# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

all return type should be list, and 
each element of the list follows the format 
(alleleIdx, [new allele1,new allele2,...])
'''

import random

__all__ = ['mutation',]

def mutation(chro, mutationRate):
    '''
    Produces a new chromosome via potential point mutation on each
    index.  If nothing changes, the original chromosome is returned.
    it is multi-points mutation.
    
    @param chro: class PyMOGEP.chromosome
    @param mutationRate, positive float: mutation mutationRate per locus
    @return: new chromosome (or self)
    '''
    genes = list(chro.genes)
    
    for geneIdx, gene in enumerate(chro.genes):
        replacements = []
        for alleleIdx, allele in enumerate(gene):
            if random.random() < mutationRate:
                if alleleIdx >= chro.headLength:
                    newAllele = random.choice(chro.terminals)
                else:
                    newAllele = random.choice(chro.symbols)
                
                # Only use this if the mutation actually did something
                if newAllele != allele:
                    replacements.append((alleleIdx, [newAllele]))
                    
        if replacements:
            genes[geneIdx] = gene.modify(replacements)
        
    return chro.newInstance(genes)

