#-*-coding:utf-8-*-
'''
Common trigonometric functions:
    - (SIN ) sin_op:         math.sin(x)
    - (COS ) cos_op:       math.cos(x)
    - (TAN ) tan_op:      math.tan(x)
    - (CSC ) csc_op:     1 / math.sin(x)
    - (SEC ) sec_op:       1 / math.cos(x)
    - (COT ) cot_op:    1 / math.tan(x)
    - (ASIN) asin_op:      math.asin(x)
    - (ACOS) acos_op:    math.acos(x)
    - (ATAN) atan_op:   math.atan(x)
    - (ACSC) acsc_op:  1 / math.asin(x)
    - (ASEC) asec_op:    1 / math.acos(x)
    - (ACOT) acot_op: 1 / math.atan(x)

Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol
import math


__all__ = 'TRIGONOMETRY_ALL', 'TRIGONOMETRY_ARITY_1'

def sin(x):
    return math.sin(x)

def cos(x):
    return math.cos(x)

def tan(x):
    return math.tan(x)

def csc(x):
    return 1./math.sin(x)

def sec(x):
    return 1./math.cos(x)

def cot(x):
    return 1./math.tan(x)

def arcsin(x):
    return math.asin(x)

def arccos(x):
    return math.acos(x)

def arctan(x):
    return math.atan(x)

def arccsc(x):
    return 1./math.asin(1)

def arcsec(x):
    return 1./math.acos(1)

def arccot(x):
    return 1./math.atan(x)


sin_op = symbol('SIN' )(sin)
cos_op = symbol('COS' )(cos)
tan_op = symbol('TAN' )(tan)
csc_op = symbol('CSC' )(csc)
sec_op = symbol('SEC' )(sec)
cot_op = symbol('COT' )(cot)
asin_op = symbol('ASIN')(arcsin)
acos_op = symbol('ACOS')(arccos)
atan_op = symbol('ATAN')(arctan)
acsc_op = symbol('ACSC')(arccsc)
asec_op = symbol('ASEC')(arcsec)
acot_op = symbol('ACOT')(arccot)

TRIGONOMETRY_ARITY_1 = sin_op, cos_op, tan_op, csc_op, \
                       sec_op, cot_op, asin_op, acos_op, \
                       atan_op, acsc_op, asec_op, \
                       acot_op
TRIGONOMETRY_ALL = TRIGONOMETRY_ARITY_1
