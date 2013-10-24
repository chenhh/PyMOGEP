#-*-coding:utf-8-*-
'''
Common exponential and logarithmic non-terminal functions and symbols:
    - (LN   ) ln_op:        math.log(x)
    - (LOG10) log10_op:     math.log10(x)
    - (^    ) power_op:     x ** y
    - (E^   ) exp_op:       e ** x
    - (10^  ) pow10_op:     10 ** x
    - (^2   ) square_op:    x ** 2
    - (^3   ) cube_op:      x ** 3
    - (Q    ) root_op:      math.sqrt(x)
    - (Q3   ) cube_root_op: x ** (1./3)
    - (^-1  ) inverse_op:   1 / x
Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol
import math

__all__ = 'POWER_ALL', 'POWER_ARITY_1', 'POWER_ARITY_2'

def ln(x):
    return math.log(x)

def log10(x):
    return math.log10(x)

def power(x, y):
    return x**y

def exp(x):
    return math.exp(x)

def power10(x):
    return 10 ** x

def square(x):
    return x*x

def cube(x):
    return x*x*x

def root(x):
    return math.sqrt(x)

def cube_root(x):
    return x ** (1./3)

def inverse(x):
    return 1./x

ln_op        = symbol('LN'   )(ln)
log10_op     = symbol('LOG10')(log10)
power_op     = symbol('^'    )(power)
exp_op       = symbol('E^'   )(exp)
pow10_op     = symbol('10^'  )(power10)
square_op    = symbol('^2'   )(square)
cube_op      = symbol('^3'   )(cube)
root_op      = symbol('Q'    )(root)
cube_root_op = symbol('Q3')(cube_root)
inverse_op   = symbol('^-1'  )(inverse)


POWER_ARITY_1 = ln_op, log10_op, root_op, exp_op, pow10_op, square_op, \
                cube_op, inverse_op
POWER_ARITY_2 = power_op,
POWER_ALL = POWER_ARITY_1 + POWER_ARITY_2
