#-*-coding:utf-8-*-
'''
Common constant functions:
    - (0) zero_op: 0
    - (1) one_op:  1
    - (P) pi_op:   math.pi
    - (E) e_op:    math.e

Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol
import math


__all__ = 'CONSTANTS_ALL', 'CONSTANTS_ARITY_0'

def zero():
    return 0

def one():
    return 1

def pi():
    return math.pi

def e():
    return math.e

zero_op = symbol('0')(zero)
one_op  = symbol('1')(one)
pi_op   = symbol('P')(pi)
e_op    = symbol('E')(e)


CONSTANTS_ARITY_0 = zero_op, one_op, pi_op, e_op
CONSTANTS_ALL = CONSTANTS_ARITY_0
