#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol
import numpy as np



@symbol('0')
def op_zero():
    return 0

@symbol('1')
def op_one():
    return 1

@symbol('pi')
def op_pi():
    return np.pi

@symbol('e')
def op_exp():
    return np.e

__all__ = [op_zero, op_one, op_pi, op_exp]