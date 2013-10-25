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

__all__ = ['inversion', 'transposeIS', 'transposeRIS', 'transposeGene']

def inversion(chro, inversionRate):
    '''
    Produces a new chromosome via headLength inversion
    @param chro: PyMOGEP.chromosome
    @param inversionRate: positive float
    @return:  chroromosome
    '''
    assert 0< inversionRate <=1.
    
    if chro.headLength < 2 or random.random() >= inversionRate: 
        return chro
    else:
        genes = list(chro.genes)
    
        geneIdx = random.choice(xrange(len(chro.genes)))
        idx1, idx2 = random.sample(xrange(chro.headLength), 2)
        if idx1 > idx2:
            idx1, idx2 = idx2, idx1
    
        # Create the new chro
        replacement = list(reversed(genes[geneIdx][idx1:idx2+1]))
        genes[geneIdx] = genes[geneIdx].modify( [idx1, replacement ])
        return chro.newInstance(genes)

def transposeIS(chro, length, transpositionISRate):
    '''
    Produces a new chromosome via IS transposition
    @param length: sequence length (typically 1, 2, or 3)
    @return: a new chromosome
    '''
    assert 0 < transpositionISRate <=1
    
    if chro.headLength < 2 or random.random() >= transpositionISRate:
        return chro
    else:
        # Pick source and target genes
        genes  = list(chro.genes)
        source = random.choice(genes)   
        target = random.choice(xrange(len(genes)))
    
        # Extract a transposition sequence. Truncate if required.
        start = random.choice(xrange(len(source)))
        end   = start + length
        end   = chro.headLength if end > chro.headLength else end
    
        # Offset into target gene: in the headLength but not the root
        offset = random.choice(xrange(1, chro.headLength))
    
        # Insert into the target gene's headLength
        replacement = source[start:end][:chro.headLength-offset] + \
                      genes[target][offset:chro.headLength-end+start]
        genes[target] = genes[target].modify([start, replacement])
        return chro.newInstance(genes)

def transposeRIS(chro, length, transpositionRISRate):
    '''
    Produces a new chro via RIS transposition
    @param length: sequence length (typically 1, 2, or 3)
    @return:       child chro
    '''
    assert transpositionRISRate
    if random.random() >= transpositionRISRate:
        return chro
    else:
        # Pick source and target genes
        genes  = list(chro.genes)
        source = random.choice(genes)
        target = random.choice(xrange(len(genes)))
        
        # Extract a transposition sequence. Truncate if required.
        # For RIS the sequence must begin with a function.
        try:
            start = random.choice(
                [idx for idx in xrange(len(source)) if callable(source[idx])]
            )
        except IndexError: # no functions!
            return chro
        
        end = start + length
        end = chro.headLength if end > chro.headLength else end
        
        # Insert into the target gene's headLength at position 0
        replacement   = source[start:end] + genes[target][:chro.headLength+start-end]
        genes[target] = genes[target].modify([start, replacement])
        return chro.newInstance(genes)

def transposeGene(chro, transpositionGeneRate):
    '''
    Produces a new chro via gene transposition
    @return: child chro
    '''
    assert transpositionGeneRate
    if len(chro.genes) < 2:
        return chro
    
    if random.random() >= transpositionGeneRate:
        return chro
    else:
        genes = list(chro.genes)
        which = random.randint(1, len(genes)-1)
        
        # Switch these genes
        genes[0], genes[which] = genes[which], genes[0]
        return chro.newInstance(genes)
