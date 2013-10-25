# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.evolution.linker import defaultLinker
from PyMOGEP.decorator import cache
import random
import itertools
from PyMOGEP.gene import Gene


class MetaChromosome(type):
    '''
    the meta class only called when first evaluate class.
    If the content is defined in class's __new__ method, 
    every time when the class is initialized, it will execute once,
    and it's very time consuming.
    
    Sets the following attributes on a chromosome class:
        - arity: maximum functional arity
        - symbols: symbols that can reside in the headLength
    Also turns caching of fitness values on for all chromosomes.
    '''
    
    def __new__(self, name, bases, dct):
        '''
        Prepares a chromosome type for use in GEP, assigning to 
        cls.symbols, cls.arity, and caching the cls._fitness.
       
        @param self: class to apply the metaclass to
        @param name: name of the class
        @param bases: base classes
        @param dct: class dict
        '''
        typ = type.__new__(self, name, bases, dct)
        typ.symbols = typ.functions + typ.terminals
        
        # Find the max arity of functions
        try:
            typ.arity = max(f.func_code.co_argcount for f in typ.functions)
        except ValueError:
            typ.arity = 0

        # Cache fitness values in attribute _ChromosomeName_cache_ (decorator)
        typ._fitnesses = cache(typ._fitnesses)
        return typ


class Chromosome(object):
    '''
    A Chromosome must provide these attributes:
        - functions: tuple of nonterminals
        - terminals: tuple of terminal symbols

    And override these functions:
        - _fitness: fitness of a given individual
        - _solved:  True if the problem is optimally solved (optional)
    '''
    __metaclass__ = MetaChromosome  # setting meta class
    _id_counter = 1   # recording current id number
    gene_type = Gene  # setting gene type

    functions = ()  # user specified
    terminals = ()  # user specified
    symbols = ()  # overridden by meta chromosome class
    headLength = tail = length = arity = 0
  
    # Unique ID of the chromosome
    chromosomeID = property(lambda self: self._id, doc='Chromosome #')
    solved = property(lambda self: self._solved(), doc='Problem solved')
    
    # multi-objective, cached in MetaChromosome
    fitnesses = property(lambda self: self._fitnesses(), 
                         doc='fitnesses')
    
    n_objectives = property(lambda self: len(self.fitnesses[0]), 
                               doc='number of objectives')
     
    @classmethod
    def randomChromosome(cls, headLength, numOfGenes=1, linker=defaultLinker):
        # class method for randomly generate genes of the chromosome
        tail = headLength * (cls.arity - 1) + 1

        newGenes = [None] * numOfGenes
        for idx in xrange(numOfGenes):
            headAlleles = [random.choice(cls.symbols)   for _ in xrange(headLength)]
            tailAlleles = [random.choice(cls.terminals) for _ in xrange(tail)]            
            newGenes[idx] = cls.gene_type(headAlleles + tailAlleles, headLength)
        return cls(newGenes, headLength, linker)

    def __init__(self, genes, headLength, linker=defaultLinker):
        '''
        @param genes: list of genes in the chro
        @param headLength: integer, length (not index) of the gene heads
        @param linker: linker function for gene evaluation
        '''
        if headLength < 0:
            raise ValueError('Head length must be >= 0')
        if not genes:
            raise ValueError('Must have at least 1 gene')

        self.genes = genes
        self.headLength = headLength
        self.linker = linker
       
        self.ParetoRank = None       # setting in evolution
        self.crowdingDistance = None # setting in evolution
        self.dominatedCount = None   # setting in evolution
        self.dominatingSet = None    # setting in evolution
     
        self._id = type(self)._id_counter
        type(self)._id_counter += 1

    def dominating(self, other):
        '''
        dominating comparison (not partial order comparison)
        @param other: other chromosome
        @return: cmp value of two chromosomes by fitness
        
        if self is dominating other then return True
        else return False
        '''
        # multi-objective (minimum is better)   
        allSame = True
        for val1, val2 in itertools.izip(self.fitnesses, other.fitnesses):
            if val1 > val2:
                return False
            if val1 != val2:
                allSame = False
        return True if not allSame else False

    def __len__(self):
        '''@return: total number of alleles in the chromosome'''
        return sum(len(g) for g in self.genes)


    def __iter__(self):
        '''@return: generator for alleles in chromosome'''
        for gene in self.genes:
            for allele in gene:
                yield allele    


    def __getitem__(self, i):
        '''
        Returns a given allele by index
        why genes[jdx][idx]
        because jdx represents the jdx-th gene
        and idx represents the idx-th allele in this gene
        @param i: allele index
        @return:  allele
        '''
        idx, jdx = divmod(i, len(self.genes))
        return self.genes[jdx][idx]


    def genesEvalRegion(self):
        '''@return: eval region of each gene (list of string)
        '''
        return [g.evalRepr for g in self.genes]
    
    def genesEvalLength(self):
        return sum(g.evalLength for g in self.genes)
    
    
    @cache
    def __repr__(self):
        '''@return: repr of chromosome alleles'''
        
        chro_str = 'ID:[%s]:%s genes:%s'%(
            self.chromosomeID,
            ''.join(repr(g) for  g in self.genes),
            ''.join("[%s]:%s "%(idx, g.reprCodingRegion()) 
                                for idx, g in enumerate(self.genes))
            )
        return chro_str


    def newInstance(self, genes):
        '''
        Creates a child chromosome
        @param genes: ordered list of GEP genes
        @return: a child chromosome of self
        '''
        return self if genes == self.genes \
            else type(self)(genes, self.headLength, self.linker)


    def eval(self, df):
        '''
        Evaluates a given GEP chromosome against some instance.  
        The terminals in the chromosome are assumed to be attributes on
        the object insgeneovided (unless they are numeric constants).
        
        @param df, pandas.DataFrame, dataset given by user
        @return: numpy.array, results of evaluating the chromosome
        '''
        return self.linker(*[g.eval(df) for g in self.genes])
    
    def _fitnesses(self):
        '''fitness function'''
        raise NotImplementedError('Must override Chromosome._fitness function')

    def _solved(self):
        '''
        It can be override by user-defined function
        @return: boolean indicating optimal solution found
        '''
        return False
