# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

C. Ferreira, "Gene Expression Programming: A New Adaptive Algorithm for 
Solving Problems.," Complex Systems, vol. 13, pp. 87-129, 2001.

prefix gene:
Xin Li, Chi Zhou, Weimin Xiao, and Peter C. Nelson. "Prefix Gene Expression 
Programming". Late Breaking Paper at the Genetic and Evolutionary Computation 
Conference(GECCO-2005), Washington, D.C., 2005.
'''
import copy
import random
from PyMOGEP.decorator import cache

class Gene(object):
    '''GEP gene class, a chromosome can contain many genes, and they 
    are linked by the linker function
    '''

    def __init__(self, alleles, headLength, RNCGenerator=None):
        '''
        @param alleles: list, symbol of function and terminal
        @param headLength: integer, head length of the gene
        @param RNCGenerator, random number generator for RNC algorithm
        
        @ivar evalLength: integer, number of alleles used after evaluating the gene
        @ivar _evalAlleles: list,  partial elements of the alleles
        @ivar evalRepr: representation of the evaluting alleles, not whole alleles   
        '''
        self.alleles = alleles
        self.headLength = headLength
        self.tailLength = len(self.alleles) - headLength
        self.RNCGenerator = RNCGenerator
        if RNCGenerator:
            self.Dc = [RNCGenerator() for _ in xrange(self.tailLength)] 
        
        self.evalLength = 0
        self._evalAlleles = None
        self._evalLength()
        self._legalForm()
        
        self.evalRepr = None
        self.evalResultArr = None
    
    def _legalForm(self):
        '''check if the head of gene is legal'''
        if self.evalLength == 0:
            self._evalLength()
        
        for allele in self.alleles[self.headLength:]:
            if callable(allele):
                errmsg = 'illegal gene, gene length:%s, head length:%s, alleles:%s'%(
                        len(self.alleles), self.headLength, self.alleles)
                raise ValueError(errmsg)
            
            
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
        

    def eval(self, df):
        '''
        Evaluates the gene against gene from argument.
        if the gene is evaluated once, the result will memory in its 
        attribute.
        we memory the evaluted the results in its attributes.
        alleles: [+, x, +, y, x ]
        => 
        1.  +
        2. x +
        3.  y x
        
        @param df, pandas.DataFrame, user specified data set
        @return, numpy.array, results of evaluating the gene.
        '''
        self._evalLength()

        jdx = self.evalLength
        for idx in reversed(xrange(jdx)):
            
            allele = self._evalAlleles[idx]
            #replace variable to corresponding array
            if isinstance(allele, str):
                if allele != '?':
                    self._evalAlleles[idx] = df[allele]
                else:
                    self._evalAlleles[idx] = self.Dc[idx - self.headLength - 1]
                    
            # if allele is a function
            if callable(allele):
                arity = allele.func_code.co_argcount
                args = self._evalAlleles[jdx - arity: jdx]
            
                # Replace the operation in eval with its return array
                self._evalAlleles[idx] = allele(*args)
                jdx -= arity

        self.evalResultArr = self._evalAlleles[0] 
        return self.evalResultArr
    
                 
    def modify(self, changes):
        '''
        if the changes are the same with alleles of this gene, 
        return this gene, not a copy.
        
        @param changes, list of (start idx, list of modify alleles)
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
            gene = type(self) (new, self.headLength, self.RNCGenerator)
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
        note: after we modify a new gene, we must update __repr__,
        else it will show old cache data, not current data.
        
        @return: string, content of gene alleles'''
        
        gene_str = ''
        for idx, allele in enumerate(self.alleles):
            try:
                # get symbol of function (added by symbol decorator)
                name = allele.symbol
            except AttributeError:
                try:
                    #get function name
                    name = allele.__name__
                except AttributeError:
                    if allele == '?':
                        name = self.Dc[idx-self.headLength-1]  
                    else: 
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



class PrefixGene(Gene):
    '''
    the eval length is the same with gene, because both the parsed trees
    have the same number of internal (function) node except 
    these functions are in different locations in the trees.
    '''
    
    def __init__(self, alleles, headLength, RNCGenerator=None ):
        super(PrefixGene, self).__init__(alleles, headLength, RNCGenerator)
        
    def eval(self, df):
        '''
        alleles: [+, x, +, y, x ]
        '''
        self._evalLength()
        
        evalStack = []
        
        arity, varCount = [], 0
        for idx, allele in enumerate(self._evalAlleles):
            if callable(allele):
                #allele is a function
                arity.append(allele.func_code.co_argcount)
                evalStack.append(allele)
            
            if isinstance(allele, str):
                #allele is a variable 
                if allele != '?':
                    allele = df[allele]
                else:
                    allele = self.Dc[idx - self.headLength - 1]
                varCount += 1
                evalStack.append(allele)
        
            while varCount == arity[-1]:
                operands = [evalStack.pop() for _ in xrange(arity)]
                operator = evalStack.pop()
                evalStack.append(operator(*operands))
                
                #use two variables to get a value
                varCount = varCount - len(operands) + 1
                arity.pop() 
                

def testGene():
    import pandas as pd
    from PyMOGEP.function.arithmetic import *
    from PyMOGEP.function.power import *
    expr = [op_root, op_multiply, op_substract, op_add, 'x', 'x', 'y', 'y' ]
    
    g = Gene([ op_substract, 'y', op_add, 'x', 'x', 'x', 'x'], 3)
    g2 = PrefixGene([ op_substract, 'y', op_add, 'x', 'x', 'x', 'x'], 3)
    print "Gene:%s, prefix: %s"%(g, g2)
    print "eval Gene:%s, prefix: %s"%(g.evalRepr, g2.evalRepr)
    print g.evalLength, g2.evalLength
  

#     x = [1,2,3]
#     y = [4,5,6]
#     df = pd.DataFrame.from_dict({"x": x, "y": y})
#     print "df:\n", df
#     print "eval:\n",g.eval(df)            

     
if __name__ == '__main__':
    testGene()
