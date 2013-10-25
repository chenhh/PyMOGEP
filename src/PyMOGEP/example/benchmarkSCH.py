#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
'''
from PyMOGEP.chromosome import Chromosome
from PyMOGEP.population import Population
from PyMOGEP.functions.mathematical.arithmetic import (add_op, subtract_op, 
    multiply_op, divide_op)
from PyMOGEP.operator.linkers import sumLinker
import random
from time import time

class Point(object):
    '''
    Schaffer's study for symbolic regression
    '''
    SAMPLE = []
    SAMPLE_SIZE = 200
    RANGE_LOW, RANGE_HIGH = -1e3, 1e3
    RANGE_SIZE = RANGE_HIGH - RANGE_LOW
    
    def __init__(self, x):
        self.x = float(x)
        
        #f1(x) = (x-2)**2
        self.y1 = (x - 2.0) * (x - 2.0)
        
        # f2(x) = x*x
        self.y2 = x * x
        
        #f3(x) = x-1
        self.y3 = x-1
        
    @staticmethod
    def populate():
        # Creates a random sample of data points
        Point.SAMPLE = []
        for _ in xrange(Point.SAMPLE_SIZE):
            x = Point.RANGE_LOW + (random.random() * Point.RANGE_SIZE)
            Point.SAMPLE.append(Point(x))

    def __repr__(self):
        return "Point:({0},[{1}, {2} {3}])".format(self.x, self.y1, self.y2, self.y3)


class SymbolicRegression(MOChromosome):

    functions = add_op, subtract_op, multiply_op, divide_op
    terminals = 'x', 
#    terminals += tuple(range(100))
    
    def _fitnesses(self):
        '''
        fitness function of the problem,
        one objective
        minimize the error1 and error2
        minize the length of gene
        '''
        error1 = 0.0
        error2 = 0.0
        error3 = 0.0
        gene1Len = float('inf')
        gene2Len = float('inf') 
        for point in Point.SAMPLE:
            try:
                guess = self(point) # Evaluation of this chromosome
#                print "guess:", guess
                error1 += abs(guess[0] - point.y1)
                error2 += abs(guess[1] - point.y2)
                error3 += abs(guess[2] - point.y3)
#                gene1Len, gene2Len = [len(gene) for gene in self.genesCodingRegion()]
#            except AttributeError: # programmer error1
#                raise
            
            except: # unviable organism
                error1 = float('inf')
                error2 = float('inf')
                error3 = float('inf')
        
        return (error1, )
    
    def _solved(self):
        '''terminal condition of this chromosome'''
        return False
    

def GEPAlgorithm(generations=10, populationSize=1000, headLength=4, numOfGenes=2):
    Point.populate()
    t1 = time()
    # Search for a solution
    p = MOPopulation(SymbolicRegression, populationSize, headLength, numOfGenes)
    p.solve(generations)
        
    #print final result
#    print "best front fitnesses:", 
#    for chro in p.bestFront:
#        print "fitness:{0}, {1}, {2}".format(chro.fitnesses, chro.genes[0].codingLength(), chro.genes[1].codingLength())
#        print "gene[0]:{0}".format(chro.genes[0].reprCodingRegion())
#        print "gene[1]:{0}".format(chro.genes[1].reprCodingRegion())
#    print "total elapsed time:{0:.3f} secs".format(time()-t1)

if __name__ == '__main__':
    
#    print fileDir
#    print packageDir
    
#    import cProfile
#    cProfile.run('GEPAlgorithm()')
    GEPAlgorithm()
