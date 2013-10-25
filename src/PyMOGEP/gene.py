# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

C. Ferreira, "Gene Expression Programming: A New Adaptive Algorithm for 
Solving Problems.," Complex Systems, vol. 13, pp. 87-129, 2001.
'''
import copy
from PyMOGEP.decorator import (memory, cache)

class Gene(object):
    '''GEP gene class, a chromosome can contain many genes, and they 
    are linked by the linker function
    '''

    def __init__(self, alleles, headLength, Dc=None):
        '''
        @param alleles: list, symbol of function and terminal
        @param headLength: integer, head length of the gene
        @param Dc: None or list, numerical constants
        
        @ivar evalLength: integer, number of alleles used after evaluating the gene
        @ivar _evalAlleles: list,  partial elements of the alleles
        @ivar evalRepr: representation of the evaluting alleles, not whole alleles   
        '''
        self.alleles = alleles
        self.headLength = headLength
        self.Dc = Dc
        
        self.evalLength = 0
        self._evalAlleles = None
        self._evalLength()
        self.evalRepr = None
        
    
    def _evalLength(self):
        '''
        the gene is parsed level order way, 
        computing number of required alleles for parsing the gene
        '''
        endIdx = 0
        for idx in xrange(self.headLength):
            if callable(self.alleles[idx]):
                endIdx += self.alleles[idx].func_code.co_argcount
            if idx == endIdx:
                break
        
        self.evalLength = endIdx + 1
        self._evalAlleles = self.alleles[:self.evalLength]
        

    @memory
    def eval(self, df, Dc):
        '''
        Evaluates the gene against gene from argument.
        if the gene is evaluated once, the result will memory in its 
        attribute.
        we memory the evaluted the results in its attributes.
        alleles: [+, x, +, y, x ]
        
        @param df, pandas.DataFrame, user specified data set
        @return, numpy.array, results of evaluating the gene.
        '''
        self._evalLength()
        self._evalLeafNodes()
        self._evalGene()
        
        idx = self.evalLength + 1
        for jdx in reversed(xrange(idx)):
            allele = self._evalAlleles[jdx]

            # if allele is a function
            if callable(allele):
                arity = allele.func_code.co_argcount
                args = self._evalAlleles[idx - arity:idx]

                #replace args to array
                arrs = [df[arg] for arg in args]
          
                # Replace the operation in eval with its return array
                self._evalAlleles[jdx] = allele(*arrs)
                idx -= arity

        return self._evalAlleles[0]
    
                 
    def derive(self, changes):
        '''
        if the changes are the same with alleles of this gene, 
        return this gene, not a copy.
        
        @param changes, list, (start idx, list of modify alleles)
        @return, gene: new gene
        '''
        new = None   # new gene
        same = True  # whether or not the evaluate region is the same
        
        for idx, alleles in changes:
            length = len(alleles)
            if self[idx: idx + length] != alleles:
                # Copy the alleles on first change
                if not new:
                    new = self[:idx] + alleles + self[idx + length:]
                else:
                    new[idx:idx + length] = alleles
                
                # Does this change the evalLength region?
                if same and idx <= self.evalLength:
                    same = False
        
        # the gene is not changed
        if not new:
            return self
        
        # the gene is changed, and it modifies the evaluate region
        if not same:
            gene = type(self) (new, self.headLength)
        # the gene is changed, but it not modifies the evaluate region
        else:
            gene = copy.copy(self)
            gene.alleles = new
            
            #update representation of the gene
            gene.evalRepr = None
            try:
                delattr(gene, "___repr___cache")
            except AttributeError:
                pass
        
        return gene
    
    @cache
    def __repr__(self):
        '''
        note: after we derive a new gene, we must update __repr__,
        else it will show old cache data, not current data.
        
        @return: string, content of gene alleles'''
        
        gene_str = ''
        for idx, allele in enumerate(self.alleles):
            # get symbol of the function (add by symbol decorator)
            try:
                name = allele.symbol
            except AttributeError:
                try:
                    name = allele.__name__
                except AttributeError:
                    name = str(allele)

            # surrounding each allele with [], 
            gene_str = "".join((gene_str, '[%s]' % name))
            if idx == self.evalLength-1:
                self.evalRepr = gene_str
        return gene_str
    
        
    def __len__(self):
        '''@return: number of alleles in the gene'''
        return len(self.alleles)
    
    
    def __iter__(self):
        '''@return: iterator over gene alleles'''
        return iter(self.alleles)
    
    def __getitem__(self, idx):
        '''@return: an individual allele from the gene'''
        return self.alleles[idx]

def testGene():
    from PyMOGEP.function.arithmetic import (
            op_add, op_substract, op_multiply, op_divide)
    g = Gene([ op_add, 'y', op_add, 'x', 'x', 'x', 'x'], 3)
    print g
    print g.evalRepr
    print g.evalLength



if __name__ == '__main__':
    testGene()
