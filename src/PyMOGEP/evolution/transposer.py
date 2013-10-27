# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2


return type of each method is list of the following format:
(alleleIdx, [new allele1,new allele2,...])
'''
import random

__all__ = ['inversion', 'transposeIS', 'transposeRIS', 'transposeGene']


def inversion(chro, inversionRate):
    '''
    Produces a new chromosome via partial head inversion
    it will not produce illegel genes.
    
    @param chro: PyMOGEP.chromosome
    @param inversionRate: positive float
    @return:  chroromosome
    '''
    assert 0< inversionRate <=1.
    
    if chro.headLength < 2 or random.random() >= inversionRate: 
        return chro
    
    genes = list(chro.genes)

    geneIdx = random.choice(xrange(len(chro.genes)))
    idx1, idx2 = random.sample(xrange(chro.headLength), 2)
    if idx1 > idx2:
        idx1, idx2 = idx2, idx1

    # Create the new chromosome
    changes = list(reversed(genes[geneIdx][idx1:idx2+1]))
    genes[geneIdx] = genes[geneIdx].modify( [[idx1, changes ]])
    return chro.newInstance(genes)


def transposeIS(chro, length, transposeISRate):
    '''
    Produces a new chromosome via IS (insertion sequence) transposition
    copy a randomly chosen seq. from idx1 with given length, 
    and overwrite the seq. idx1 from idx2  of the same chromosome
    
    @param length: sequence length (typically 1, 2, or 3)
    @return: a new chromosome
    '''
    assert 0 < transposeISRate <=1.
    
    if chro.headLength < 2 or random.random() >= transposeISRate:
        return chro
   
    # Pick srcGene and tgtGene genes
    genes  = list(chro.genes)
    srcGene = random.choice(genes)   
    tgtGene = random.choice(xrange(len(genes)))

    # Extract a transposition sequence.
    idx1 = random.choice(xrange(len(srcGene)))
    idx2   = idx1 + length
    idx2   = chro.headLength if idx2 > chro.headLength else idx2

    # Offset into tgtGene gene: in the headLength but not the root
    offset = random.choice(xrange(1, chro.headLength))

    # Insert into the tgtGene gene's headLength
    changes = srcGene[idx1:idx2][:chro.headLength-offset] + \
                  genes[tgtGene][offset:chro.headLength-idx2+idx1]
    genes[tgtGene] = genes[tgtGene].modify([[idx1, changes]])
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
        genes[target] = genes[target].modify([[start, replacement]])
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
