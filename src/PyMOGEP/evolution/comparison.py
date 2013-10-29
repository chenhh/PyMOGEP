# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''
__all__ = ['partialOrder',]

def partialOrder(chro1, chro2):
    '''
    ordering by the rank and the crowding distance 
    @param chro1, PyMOGEP.chromosome
    @param chro2, PyMOGEP.chromosome 
    '''
    #same instance
    if chro1 is chro2:
        return 0
     
    #prefer better rank
    if chro1.ParetoRank < chro2.ParetoRank:
        return 1
    #larger crowing distance is better
    elif (chro1.ParetoRank == chro2.ParetoRank and 
        chro1.crowdingDistance > chro2.crowdingDistance):
        return 1
    elif  (chro1.ParetoRank == chro2.ParetoRank and 
        chro1.crowdingDistance == chro2.crowdingDistance):
        return 0
    else:
        return -1