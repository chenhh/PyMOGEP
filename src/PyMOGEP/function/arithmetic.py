#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

from PyMOGEP.decorator import symbol


@symbol('+')
def add(x, y):
    return x + y

@symbol('-')
def substract(x, y):
    return x - y

@symbol('*')
def multiply(x, y):
    return x * y

@symbol('/')
def divide(x, y):
    return float(x) / y

@symbol('%')
def modulus(x, y):
    return x % y
 
__all__ = add, substract, multiply, divide, modulus
