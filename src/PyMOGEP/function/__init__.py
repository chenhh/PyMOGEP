#-*-coding:utf-8-*-

'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.function.arithmetic import *
from PyMOGEP.function.comparison import *
from PyMOGEP.function.constants import *
from PyMOGEP.function.hyperbolic import *
from PyMOGEP.function.power import *
from PyMOGEP.function.rounding import *
from PyMOGEP.function.trigonometry import *


__all__ = 'MATH_ALL', 'MATH_ARITY_0', 'MATH_ARITY_1', 'MATH_ARITY_2'


MATH_ARITY_0 = CONSTANTS_ARITY_0
MATH_ARITY_1 = HYPERBOLIC_ARITY_1 + POWER_ARITY_1 + ROUNDING_ARITY_1 + \
               TRIGONOMETRY_ARITY_1
MATH_ARITY_2 = ARITHMETIC_ARITY_2 + COMPARISON_ARITY_2 + POWER_ARITY_2
MATH_ALL = MATH_ARITY_0 + MATH_ARITY_1 + MATH_ARITY_2
