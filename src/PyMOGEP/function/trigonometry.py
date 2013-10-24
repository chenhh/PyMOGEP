#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol
import numpy as np



@symbol('SIN')
def op_sin(x):
    return np.sin(x)

@symbol('COS')
def op_cos(x):
    return np.cos(x)

@symbol('TAN')
def op_tan(x):
    return np.tan(x)

@symbol('CSC')
def op_csc(x):
    return 1./np.sin(x)

@symbol('SEC')
def op_sec(x):
    return 1./np.cos(x)

@symbol('COT')
def op_cot(x):
    return 1./np.tan(x)

@symbol('ASIN')
def op_arcsin(x):
    return np.arcsin(x)

@symbol('ACOS')
def op_arccos(x):
    return np.arccos(x)

@symbol('ATAN')
def op_arctan(x):
    return np.arctan(x)

@symbol('ACSC')
def op_arccsc(x):
    return 1./np.arcsin(x)

@symbol('ASEC')
def op_arcsec(x):
    return 1./np.arccos(x)

@symbol('ACOT')
def op_arccot(x):
    return 1./np.arctan(x)

@symbol('SINH')
def op_sinh(x):
    return np.sinh(x)

@symbol('COSH')
def op_cosh(x):
    return np.cosh(x)

@symbol('TANH')
def op_tanh(x):
    return np.tanh(x)

@symbol('CSCH')
def op_csch(x):
    return 1./np.sin(x)

@symbol('SECH')
def op_sech(x):
    return 1./np.cosh(x)

@symbol('COTH')
def op_coth(x):
    return 1./np.tanh(x)

__all__ = [op_sin, op_cos, op_tan, op_csc, op_sec, op_cot,
           op_arcsin, op_arccos, op_arctan, op_arccsc, op_arcsec, op_arccot,
           op_sinh, op_cosh, op_tanh, op_csch, op_sech, op_coth]