#-*-coding:utf-8-*-
'''
Provides basic comparison non-terminals.  If each is true, it returns
the first value given.  If false, it returns the second.

Common comparison non-terminal functions:
    - (=) equal_op:            i if i == j else j
    - (U) unequal_op:          i if i != j else j
    - (<) less_op:             i if i < j else j
    - (>) greater_op:          i if i > j else j
    - (L) less_or_equal_op:    i if i <= j else j
    - (G) greater_or_equal_op: i if i >= j else j

Note:
使用lambda為算式，會造成物件無法pickle,而無法做multi processing
所以改用def function
'''

from PyMOGEP.functions import symbol

__all__ = 'COMPARISON_ALL', 'COMPARISON_ARITY_2'

#def equal(x, y):
#    return x if x == y else y
#
#def unequal(x, y):
#    return x if x != y else y
#
#def less(x, y):
#    return x if x < y else y
#
#def greater(x, y):
#    return x if x > y else y
#    
#def less_or_equal(x, y):
#    return x if x <= y else y
#
#def greater_or_equal(x, y):
#    return x if x>=y else y

def equal(x, y):
    return x==y
def unequal(x, y):
    return x != y 

def less(x, y):
    return x < y

def greater(x, y):
    return x > y
    
def less_or_equal(x, y):
    return x <= y

def greater_or_equal(x, y):
    return x>=y

equal_op            = symbol('=')(equal)
unequal_op          = symbol('U')(unequal)
less_op             = symbol('<')(less)
greater_op          = symbol('>')(greater)
less_or_equal_op    = symbol('L')(less_or_equal)
greater_or_equal_op = symbol('G')(greater_or_equal)


COMPARISON_ARITY_2 = equal_op, unequal_op, less_op, greater_op, \
                     less_or_equal_op, greater_or_equal_op
COMPARISON_ALL = COMPARISON_ARITY_2
