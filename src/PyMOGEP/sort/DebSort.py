# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective algorithm: NSGA-II

K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, "A fast and elitist multiobjective 
genetic algorithm: NSGA-II," Evolutionary Computation, IEEE Transactions on, 
vol. 6, pp. 182-197, 2002.
'''

def DebNonDominatedSort(population):
        '''
        @param population, PyMoGEP.population
        @ivar chromosome.dominatedCount: 
               number of the other chromosome dominating this one
        @ivar chromosome.dominatingSet: 
            the chromosomes dominated by this one
        @return all Pareto fronts (list)
        '''
        #1st front 
        ParetoFronts=[[]]
        
        #finding out the 1st Pareto front
        for idx, chro1 in enumerate(population):
            chro1.dominatedCount = 0
            chro1.dominatingSet = []
            
            for jdx, chro2 in enumerate(population):
                if idx != jdx:
                    #chro1 dominating chro2
                    if chro1.dominating(chro2):
                        chro1.dominatingSet.append(chro2)
                        
                    #chro1 dominated by chro2
                    elif chro2.dominating(chro1):
                        chro1.dominatedCount += 1
                        
            if chro1.dominatedCount == 0:
                chro1.ParetoRank = 1
                ParetoFronts[0].append(chro1)
        
        #searching other Pareto fronts
        frontIdx = 0
        while len( ParetoFronts[frontIdx] ) > 0:              
            nextFront=[]
            for chro1 in ParetoFronts[frontIdx]:
                for chro2 in chro1.dominatingSet:
                    chro2.dominatedCount -= 1                  
                    if chro2.dominatedCount == 0:
                        chro2.ParetoRank = (frontIdx + 2)
                        nextFront.append(chro2)
            frontIdx += 1
            ParetoFronts.append(nextFront)
        
        #remove last empty front
        assert len(ParetoFronts)>1 and len(ParetoFronts[-1]) == 0
        ParetoFronts.pop()

        return ParetoFronts
    
