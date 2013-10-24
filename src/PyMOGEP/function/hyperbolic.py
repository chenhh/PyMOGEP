#-*-coding:utf-8-*-
'''
Common hyperbolic functions:
    - (SINH) sinh_op:      math.sinh(x)
    - (COSH) cosh_op:    math.cosh(x)
    - (TANH) tanh_op:   math.tanh(x)
    - (CSCH) csch_op:  1 / math.sinh(x)
    - (SECH) sech_op:    1 / math.cosh(x)
    - (COTH) coth_op: 1 / math.tanh(x)
Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol
import math


__all__ = 'HYPERBOLIC_ALL', 'HYPERBOLIC_ARITY_1'

def sinh(x):
    return math.sinh(x)

def cosh(x):
    return math.cosh(x)

def tanh(x):
    return math.tanh(x)

def csch(x):
    return 1./math.sin(x)

def sech(x):
    return 1./math.cosh(x)

def coth(x):
    return 1./math.tanh(x)

sinh_op = symbol('SINH')(sinh)
cosh_op = symbol('COSH')(cosh)
tanh_op = symbol('TANH')(tanh)
csch_op = symbol('CSCH')(csch)
sech_op = symbol('SECH')(sech)
coth_op = symbol('COTH')(coth)


HYPERBOLIC_ARITY_1 = sinh_op, cosh_op, tanh_op, csch_op, \
                     sech_op, coth_op
HYPERBOLIC_ALL = HYPERBOLIC_ARITY_1
