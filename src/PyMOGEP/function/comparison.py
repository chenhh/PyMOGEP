#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol


@symbol('==')
def op_equal(x, y):
    return x==y

@symbol('!=')
def op_unequal(x, y):
    return x != y 

@symbol('<')
def op_less(x, y):
    return x < y

@symbol('>')
def op_greater(x, y):
    return x > y

@symbol('<=')    
def op_less_or_equal(x, y):
    return x <= y

@symbol('>=')
def op_greater_or_equal(x, y):
    return x>=y


__all__ = [op_equal, op_unequal, op_less, op_greater, op_less_or_equal, 
           op_greater_or_equal]
