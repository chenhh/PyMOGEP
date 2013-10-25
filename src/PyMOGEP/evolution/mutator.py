# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

import random


def mutation(chro, rate):
    '''
    Produces a new chro via potential point mutation on each
    locus.  If nothing changes, the original chro is returned.
    
    @param chro: class PyMOGEP.chro.Chromosome
    @param rate: mutation rate per locus
    @return:     child chro (or self)
    '''
    genes = list(chro.genes)
    
    # Traverse the chro gene by gene, then locus by locus
    for geneIdx, gene in enumerate(chro.genes):
        replacements = []
        for i, allele in enumerate(gene):
            # Do we mutate this locus?
            if random.random() < rate:
                if i >= chro.headLength:
                    newAllele = random.choice(chro.terminals)
                else:
                    newAllele = random.choice(chro.symbols)
                
                # Only use this if the mutation actually did something
                if newAllele != allele:
                    replacements.append((i, [newAllele]))

        # If we have actual replacements to make, do them
        if replacements:
            genes[geneIdx] = gene.derive(replacements)
        
    # Create a child of this chro
    return chro._child(genes)
