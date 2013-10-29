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
    idx1, idx2 = random.sample(xrange(chro.headLength+1), 2)
    if idx1 > idx2:
        idx1, idx2 = idx2, idx1

    changes = list(reversed(genes[geneIdx][idx1:idx2]))
    genes[geneIdx] = genes[geneIdx].modify( [[idx1, changes ]])
    return chro.newInstance(genes)


def transposeIS(chro, length, transposeISRate):
    '''
    Produces a new chromosome via IS (insertion sequence) transposition
    copy a randomly chosen seq. from idx1 with given length, 
    and overwrite the seq. from idx1 to idx2  of the same chromosome
    
    @param length: sequence length (typically 1, 2, or 3)
    @return: a new chromosome
    '''
    assert 0 < transposeISRate <=1.
    
    if chro.headLength < 2 or random.random() >= transposeISRate:
        return chro
   
    # Pick srcGeneIdx and tgtGeneIdx genes
    genes  = list(chro.genes)
    geneLength = len(genes)
    srcGeneIdx = random.choice(xrange(geneLength))   
    tgtGeneIdx = random.choice(xrange(geneLength))

    # Extract a transposition sequence.
    idx1 = random.choice(xrange(geneLength))
    idx2   = idx1 + length
    idx2   = chro.headLength if idx2 > chro.headLength else idx2

    #start idx of the tgtGeneIdx gene, but not the root
    tgtIdx1 = random.choice(xrange(1, chro.headLength))

    # Insert into the tgtGeneIdx gene's headLength
    changes = genes[srcGeneIdx][idx1:idx2][:chro.headLength-tgtIdx1] + \
                genes[tgtGeneIdx][tgtIdx1:chro.headLength-idx2+idx1]
    genes[tgtGeneIdx] = genes[tgtGeneIdx].modify([[idx1, changes],])
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
    
    genes  = list(chro.genes)
    geneLength = len(genes)
    srcGeneIdx = random.choice(xrange(geneLength))
    tgtGeneIdx = random.choice(xrange(geneLength))
    
    # Extract a transposition sequence. Truncate if required.
    # For RIS, the sequence must begin with a function.
    try:
        idx1 = random.choice(
            [idx for idx in xrange(geneLength) 
             if callable(genes[srcGeneIdx][idx])]
        )
    except IndexError: 
        # no functions after the idx1 of the source gene
        return chro
    
    idx2 = idx1 + length
    idx2 = chro.headLength if idx2 > chro.headLength else idx2
    
    # Insert into the tgtGeneIdx gene's headLength at position 0
    changes = genes[srcGeneIdx][idx1:idx2] + genes[tgtGeneIdx].alleles
    changes = changes[:chro.headLength] + genes[tgtGeneIdx][chro.headLength:]
    
    genes[tgtGeneIdx] = genes[tgtGeneIdx].modify([[0, changes],])
    return chro.newInstance(genes)


def transposeGene(chro, transpositionGeneRate):
    '''
    Produces a new chromosome via gene transposition
    @return: child chromosome
    '''
    assert transpositionGeneRate
    if len(chro.genes) < 2:
        return chro
    
    if random.random() >= transpositionGeneRate:
        return chro
    else:
        genes = list(chro.genes)
        idx = random.randint(1, len(genes)-1)
        
        # Switch these genes
        genes[0], genes[idx] = genes[idx], genes[0]
        return chro.newInstance(genes)
