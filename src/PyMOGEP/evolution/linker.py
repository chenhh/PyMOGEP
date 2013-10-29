# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

For combining multi-genetic chromosome.
'''


__all__ = ['defaultLinker', 'sumLinker', 'andLinker', 'orLinker']


def defaultLinker(*args):
    '''@return: either a single value or a tuple, depending on context'''
    if len(args) == 1:
        return args[0]
    return args


def sumLinker(*args):
    '''@return: the sum of all sub-ETs'''
    return sum(args)


def andLinker(*args):
    '''@return: the AND of all given args'''
    return all(args)


def orLinker(*args):
    '''@return: the OR of all given args'''
    return any(args)
