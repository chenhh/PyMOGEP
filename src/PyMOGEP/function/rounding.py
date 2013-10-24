#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol
import numpy as np



@symbol('FLOOR')
def op_floor(x):
    return np.floor(x)

@symbol('CEIL')
def op_ceil(x):
    return np.ceil(x)

@symbol('ROUND')
def op_round(x):
    return np.round(x)

@symbol('ABS')
def op_abs(x):
    return abs(x)

__all__ = [op_floor, op_ceil, op_round, op_abs]