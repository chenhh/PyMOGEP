#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw


K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, "A fast and elitist multiobjective 
genetic algorithm: NSGA-II," Evolutionary Computation, IEEE Transactions on, 
vol. 6, pp. 182-197, 2002.

number of variable: 1
variables bounds: [-10**3, 10**3]
objective functions:
    1. f(x) = x**2
    2. g(x) = (x-2)**2

optimal solution: x in [0, 2]
'''
from PyMOGEP.chromosome import Chromosome
from PyMOGEP.population import Population
from PyMOGEP.function.arithmetic import *
from PyMOGEP.evolution.linker import *
import random
import pandas as pd
import numpy as np
from time import time


def Dataset(n_data=1000):
    func0 = lambda x: x**2
    func1 = lambda x: (x-2)**2
    lower, upper = -10**3, 10**3
    x = (upper-lower) * np.random.random((n_data)) + lower
    f1 =func0(x)
    f2 =func1(x)
    return pd.DataFrame.from_dict({"x": x, 
                                   "f1": f1,
                                   "f2": f2}) 


class SymbolicRegression(Chromosome):

    functions = op_add, op_multiply, op_substract, op_divide
    terminals = 'x', 

    def _fitnesses(self):
        '''fitness function'''
        error1 = 0.0
        error2 = 0.0
        
        # Evaluation of this chromosome
        guess = self.eval(Population.df)

        error1 += abs(guess[0] - Population.df['f1'])
        error2 += abs(guess[1] - Population.df['f2'])
        print "e1:", error1
        print "e2:", error2
        return (error1, error2)
    
    def _solved(self):
        '''termination condition'''
        return False
    

def GEPAlgorithm(generations=10, popSize=1000, 
                 headLength=4, n_genes=2):
    df = Dataset()
    print df
    t0 = time()

    Population.df = df
    p = Population(SymbolicRegression, popSize, headLength, n_genes)
    p.solve(generations)
        
    #print final result
    print "best fitnesses:", 
    for chro in p.bestFront:
        print "fitness:%s, %s, %s"%(chro.fitnesses, 
                                    chro.genes[0].evalLength, 
                                    chro.genes[1].evalLength)
        print "gene[0]:%s"%(chro.genes[0].evalRepr)
        print "gene[1]:%s"%(chro.genes[1].evalRepr)
    print "5.3f secs"%(time()-t0)

if __name__ == '__main__':
    
#    import cProfile
#    cProfile.run('GEPAlgorithm()')
    GEPAlgorithm()
