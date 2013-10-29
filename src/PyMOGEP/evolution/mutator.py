# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

return type of each method is list of the following format:
(alleleIdx, [new allele1,new allele2,...])
'''

import random

__all__ = ['mutation',]

def mutation(chro, mutationRate):
    '''
    Produces a new chromosome via multi-point mutation on each
    index.  If nothing changes, the original chromosome is returned.
    
    @param chro: PyMOGEP.chromosome
    @param mutationRate, positive float
    @return: new chromosome (or self)
    '''
    assert 0. < mutationRate <=1.
    genes = list(chro.genes)
    
    for geneIdx, gene in enumerate(chro.genes):
        changes = []
        for alleleIdx, allele in enumerate(gene):
            if random.random() < mutationRate:
                if alleleIdx >= chro.headLength:
                    newAllele = random.choice(chro.terminals)
                else:
                    newAllele = random.choice(chro.symbols)
                
                # Only use this if the mutation actually did something
                if newAllele != allele:
                    changes.append((alleleIdx, [newAllele,]))
                    
        if changes:
            genes[geneIdx] = gene.modify(changes)
        
    return chro.newInstance(genes)

