# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''
import random


def invert(chro, inversionRate):
    '''
    Produces a new chro via headLength inversion
    @param inversionRate
    @return: child chro
    '''
    assert inversionRate
    
    if chro.headLength < 2 or random.random() >= inversionRate: 
        return chro
    else:
        genes = list(chro.genes)
    
        #idx: 那一個gene要做invert
        #start, stop: gene head invert的起始點與終點
        idx = random.choice(xrange(len(chro.genes)))
        start, stop = random.sample(xrange(chro.headLength), 2)
        if start > stop:
            start, stop = stop, start
    
        # Create the new chro
        replacement = list(reversed(genes[idx][start:stop+1]))
        genes[idx] = genes[idx].derive([(start, replacement)])
        return chro._child(genes)

def transposeIS(chro, length, transpositionISRate):
    '''
    Produces a new chro via IS transposition
    @param length: sequence length (typically 1, 2, or 3)
    @return:       child chro
    '''
    # Since IS does not transpose to the root, it has no purpose
    # if the headLength length is less than 2.
    assert transpositionISRate
    
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
        genes[target] = genes[target].derive([(offset, replacement)])
        return chro._child(genes)

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
        genes[target] = genes[target].derive([(0, replacement)])
        return chro._child(genes)

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
        return chro._child(genes)
