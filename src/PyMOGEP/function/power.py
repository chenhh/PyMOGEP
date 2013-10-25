#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol
import numpy as np

__all__ = ['op_ln', 'op_log10', 'op_power', 'op_exp', 'op_power10', 'op_square',
           'op_cube', 'op_root', 'op_cube_root', 'op_inverse']

@symbol('LN')
def op_ln(x):
    return np.log(x)

@symbol('LOG10')
def op_log10(x):
    return np.log10(x)

@symbol('^')
def op_power(x, y):
    return x**y

@symbol('E^')
def op_exp(x):
    return np.exp(x)

@symbol('10^')
def op_power10(x):
    return 10 ** x

@symbol('^2')
def op_square(x):
    return x*x

@symbol('^3')
def op_cube(x):
    return x*x*x

@symbol('Q')
def op_root(x):
    return np.sqrt(x)

@symbol('Q3')
def op_cube_root(x):
    return x ** (1./3)

@symbol('^-1')
def op_inverse(x):
    return 1./x

