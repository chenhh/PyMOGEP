# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

C. Ferreira, "Gene Expression Programming: A New Adaptive Algorithm for 
Solving Problems.," Complex Systems, vol. 13, pp. 87-129, 2001.
'''
import numpy as np
import itertools
from operator import itemgetter
from PyMOGEP.decorator import memory

class Gene(object):
    '''GEP gene class, a chromosome can contain many genes, and they 
    are linked by the linker function
    '''

    def __init__(self, alleles, headLength):
        '''
        @param alleles: numpy.array, symbol of function and terminal
        @param headLength: integer, head length of the gene
        
        @ivar evalLength: integer, number of alleles used after evaluating the gene
        @ivar _evalAlleles: numpy.array,  partial elements of the alleles
        @ivar _leafNodes: numpy.array of list, the list contains the information 
                          (variable, [location])
        
        '''
        self.alleles = alleles
        self.headLength = headLength
    
        self.evalLength = -1
        self._evalAlleles = None
        self._leafAlleles = None
        self._parseGene()
    
    def _evalLength(self):
        '''
        the gene is parsed level order way, 
        computing number of required alleles for parsing the gene
        '''
        idx, pos = 0, 1
        while pos:
            next_pos = 0
            for _ in xrange(pos):
                # computing the arity if the allele is a function
                if callable(self.alleles[idx]):
                    next_pos += self.alleles[idx].func_code.co_argcount
                idx += 1
            pos = next_pos
        self.evalLength = idx - 1
        self._evalAlleles = self.alleles[:idx]
        
    def _aggregrateLeafAlleles(self):
        # the variable in the alleles is string 
        # self._evaluation=[+, +, x, x, x, y ]
        # leafs = [(2,x), (3, x), (4, x), (5, y)]
        # aggregate the leaf nodes, (variable, [locations])
        # terminals = [('x', [2, 3]), (('y', [5])),...]
        
        leafs = np.array([(idx, content) for idx, content in 
                          enumerate(self._evalAlleles) 
                          if isinstance(content, str)])
        
        self._leafAlleles = [ (key, [idx[0] for idx in value]) for 
                             key, value in itertools.groupby(
                            sorted(leafs, key=itemgetter(1)), key=itemgetter(1))
        ]
        
    def _evalGene(self, df):
        '''
        replace the variables in self._evalAlleles to the values
        defined by the user.
        @param dataframe, pandas.DataFrame, data given by users for evaluating
                                    and it should contain one row 
        '''
        
      
    
    @memory
    def parse(self, gene):
        '''
        E.g.:
        gene = MOKarvaGene()   #call __init__
        gene()                 #call __call__
        
        Evaluates the Karva gene against gene from argument.
        if the gene is evaluated, the result will memoize in its 
        attribute.
          
        @param gene: gene of a Karva gene
        @return: result of evaluating the gene.
        '''
        self._prepareEvalAttributes(gene)
        
        # Evaluate the gene against gene in reverse
        index = self.codingLength + 1
        for i in reversed(xrange(index)):
            allele = self.alleles[i]

            # if allele is a function
            if callable(allele):
                arity = allele.func_code.co_argcount
                args = self._evaluation[index - arity:index]

                # Replace the operation in eval with its return val
                self._evaluation[i] = allele(*args)
                index -= arity

        return self._evaluation[0]
    
                 
    def derive(self, changes):
        '''
        Derives a gene from self.
        type1: (修改的位置, 修改後的符號)
        e.g. replacing allele 0 with 'x' and allele 3 with a 
        function add via point mutation would look like:
            gene.derive([(0, ['x']), (3, [add]))
        
        type2: (修改的起始點, 連續修改後的符號)
        e.g. replacing a block of alleles in crossover would like like:
            gene.derive([(5, [add, 'x', 'y'])])

        #如果changes後的部份與原來相同，則傳回原物件(not a copy)
        @param changes: sequence of (index, alleles) tuples
        @return: new KarvaGene
        '''
        new = None  # new gene
        same = True  # whether or not the codingLength region is the same
        
        for index, alleles in changes:
            length = len(alleles)
            if self[index:index + length] != alleles:
                # Copy the alleles on first change
                if not new:
                    new = self[:index] + alleles + self[index + length:]
                else:
                    new[index:index + length] = alleles
                
                # Does this change the codingLength region?
                if same and index <= self.codingLength:
                    same = False
        
        # 完全沒有修改gene
        if not new:
            return self
        
        # gene被修改,且有修改到coding region
        if not same:
            gene = type(self) (new, self.headLength)
        # gene被修改,但未修改到coding region
        else:
            gene = copy(self)
            gene.alleles = new
            gene.codingRepr = None
            try:
                delattr(gene, "___repr___cache")
            except AttributeError:
                pass
        
        return gene
    
    def __repr__():
        pass
    
        
    def __len__(self):
        '''@return: number of alleles in the gene'''
        return self.alleles.size
        
    def __iter__(self):
        '''@return: iterator over gene alleles'''
        return iter(self.alleles)
    
    def __getitem__(self, idx):
        '''@return: an individual allele from the gene'''
        return self.alleles[idx]

if __name__ == '__main__':
    pass