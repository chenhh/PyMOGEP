#-*-coding:utf-8-*-
'''
Common rounding functions:
    - (FLOOR) floor_op: math.floor(x)
    - (CEIL ) ceil_op:  math.ceil(x)
    - (ROUND) round_op: round(x)
    - (ABS  ) abs_op:   abs(x)
Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol
import math


__all__ = 'ROUNDING_ALL', 'ROUNDING_ARITY_1'

def floor(x):
    return math.floor(x)

def ceil(x):
    return math.ceil(x)

def roundop(x):
    return round(x)

def absop(x):
    return abs(x)


floor_op = symbol('FLOOR')(floor)
ceil_op  = symbol('CEIL' )(ceil)
round_op = symbol('ROUND')(roundop)
abs_op   = symbol('ABS'  )(absop)


ROUNDING_ARITY_1 = floor_op, ceil_op, round_op, abs_op
ROUNDING_ALL = ROUNDING_ARITY_1
