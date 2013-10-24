'''
Created on 2013/10/24

@author: chenhh
'''

from PyMOGEP.function.arithmetic import (op_add, op_substract)
from PyMOGEP.function.constants import (op_pi, op_exp)

def requiredLen(alleles):
    endPos = 0
    for idx in xrange(len(alleles)):
        if callable(alleles[idx]):
            endPos += alleles[idx].func_code.co_argcount
        if idx == endPos:
            break
    print endPos

if __name__ == '__main__':
    alleles = ['x', op_add,  op_pi,  op_add, 'y', 'y','y','y','y','y','y','y']
    requiredLen(alleles)