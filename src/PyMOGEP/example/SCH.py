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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time


def Dataset(n_data=1000):
    func1 = lambda x: x**2
    func2 = lambda x: (x-1)**2
  
    lower, upper = -10**3, 10**3
    x = (upper-lower) * np.random.random((n_data)) + lower
    f1 = func1(x)
    f2 = func2(x)
    return pd.DataFrame.from_dict({"x": x, "f1": f1, "f2": f2}) 


class SymbolicRegression(Chromosome):

    functions = op_add, op_multiply, op_substract
    terminals = 'x', 

    def _fitnesses(self):
        '''fitness function'''

        # Evaluation of this chromosome
        guess = self.eval(Population.df)

        try:
            error = ( np.sum(np.abs(guess[0] - Population.df['f1'])),
                      np.sum(np.abs(guess[1] - Population.df['f2']))
                      )
            return error
        
        except Exception:
            return (float('inf'), float('inf'))
        
    def _solved(self):
        '''termination condition'''
        return False
    

def GEPAlgorithm(generations=10, popSize=500, 
                 headLength=4, n_genes=2):
    df = Dataset(100)
    t0 = time()
 
    Population.df = df
    p = Population(SymbolicRegression, popSize, headLength, n_genes,
                   n_elites=2, RNCGenerator=np.random.randn, verbose=False)
    p.solve(generations)
         
    #print final result
    print "best Pareto front:", 
    for chro in p.bestFront:
        print "fitness:%s"%(chro.fitnesses,)
        print chro
    print "%.3f secs"%(time()-t0)

if __name__ == '__main__':
#    import cProfile
#    cProfile.run('GEPAlgorithm()')  
    GEPAlgorithm()
