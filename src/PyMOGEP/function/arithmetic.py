#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol

__all__ = ['op_add', 'op_substract', 'op_multiply', 'op_divide', 'op_modulus']

@symbol('+')
def op_add(x, y):
    return x + y

@symbol('-')
def op_substract(x, y):
    return x - y

@symbol('*')
def op_multiply(x, y):
    return x * y

@symbol('/')
def op_divide(x, y):
    return float(x) / y

@symbol('%')
def op_modulus(x, y):
    return x % y
 