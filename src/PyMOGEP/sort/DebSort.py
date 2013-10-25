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
        @param population
        chromosome.dominatedCount: 
            the number of chromosomes dominating this chromsome
        chromosome.dominatingSet: 
            the chromosomes dominated by this chromosome
        @return all Pareto fronts (list)
        '''
        #1st front must exist
        ParetoFronts=[[]]
        
        #finding out the 1st Pareto front
        for idx, chro1 in enumerate(population):
#            print "outer non-dominated chro id:", chro1.chromosomeID, id(chro1), id(chro1.genes[0]), id(chro1.genes[1])
            chro1.dominatedCount = 0
            chro1.dominatingSet = []
            
            for jdx, chro2 in enumerate(population):
#                print "inner non-dominated chro id:", chro2.chromosomeID, id(chro2), id(chro2.genes[0]), id(chro2.genes[1])
#                print "chro1: fitness", chro1.fitnesses
#                print "chro2: fitness", chro2.fitnesses
                if idx != jdx:
                    #chro1 dominating chro2
                    if chro1.fitnessesDominating(chro2):
                        chro1.dominatingSet.append(chro2)
                        
                    #chro1 dominated by chro2
                    elif chro2.fitnessesDominating(chro1):
                        chro1.dominatedCount += 1
#                print "inner after dominating"
            if chro1.dominatedCount == 0:
                chro1.ParetoRank = 1
                ParetoFronts[0].append(chro1)
        
#        print "non-dominated best frontier"
        #searching other Pareto fronts
        ParetoIdx = 0
        while len( ParetoFronts[ParetoIdx] ) > 0:              
            nextFront=[]
            for chro1 in ParetoFronts[ParetoIdx]:
                for chro2 in chro1.dominatingSet:
                    chro2.dominatedCount -= 1                  
                    if chro2.dominatedCount == 0:
                        chro2.ParetoRank = (ParetoIdx + 2)
                        nextFront.append(chro2)
            ParetoIdx += 1
#            print "non-dominated Pareto idx:{0}".format(ParetoIdx)
            ParetoFronts.append(nextFront)
        
        #remove last empty front
        assert len(ParetoFronts)>1 and len(ParetoFronts[-1]) == 0
        ParetoFronts.pop()

        return ParetoFronts