'''
Created on 2013/10/24

@author: chenhh
'''

from PyMOGEP.function.arithmetic import (op_add, op_substract)
from PyMOGEP.function.constants import (op_pi, op_exp)
from time import time
import numpy as np
import numexpr as ne

def requiredLen(alleles):
    endPos = 0
    for idx in xrange(len(alleles)):
        if callable(alleles[idx]):
            endPos += alleles[idx].func_code.co_argcount
        if idx == endPos:
            break
    print endPos


def group(alleles):
    leafs = {}
    for idx, allele in enumerate(alleles):
        if isinstance(allele, str):
            if allele in leafs.keys():
                leafs[allele].append(idx)
            else:
                leafs[allele] = [idx,]
    
    arr = [(key, np.asarray(val)) for key, val in leafs.items()]
    print arr
    

def testExpression():
    
    n = 1000
    a = np.random.randn(n)
    b = np.random.randn(n)
    
    t0 = time()
    c = a+b-b**2
    print "np, %.3f secs"%(time()-t0)
    
    t0= time()
    ne.evaluate('a+b-b**2')
    print "ne, %.3f secs"%(time()-t0)
    
if __name__ == '__main__':
#     alleles = ['x', op_add,  op_pi,  op_add, 'y', 'y','y','y','y','y','y','y']
#     requiredLen(alleles)
#     group(alleles)
#     print np.asarray(v for v in xrange(10))
    testExpression()