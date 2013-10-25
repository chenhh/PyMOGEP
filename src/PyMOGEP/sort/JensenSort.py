# -*- coding: utf-8 -*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2

multi-objective algorithm: NSGA-II

M. T. Jensen, "Reducing the run-time complexity of multiobjective EAs: 
The NSGA-II and other algorithms," Evolutionary Computation, IEEE Transactions on, 
vol. 7, pp. 503-515, 2003.
'''

import bisect

def twoObjectivCmpFun(chro1, chro2):
    '''
    if chromosome1.fitnesses dominating chromosome2.fitnesses, then return 1
    elif chromosome1.fitnesses == chromosome2.fitnesses, then return 0
    else return -1
    '''
    assert chro1.numOfObjectives == chro2.numOfObjectives == 2
    
    if (chro1.fitnesses[0] < chro2.fitnesses[0]) or (
     chro1.fitnesses[0] == chro2.fitnesses[0] and 
     chro1.fitnesses[1] < chro2.fitnesses[1]):
        return 1
    elif chro1.fitnesses == chro2.fitnesses:
        return 0
    else:
        return -1

def twoObjectivesSweepAlgorithm(population):
    '''time complexity O(NlogN)
    M: num. of objective
    N: num. of points
    '''
    #decreasing (dominating), O(NlogN)
    population.sort(cmp=twoObjectivCmpFun, reverse=True)
    ParetoFronts =[[population[0]]]
    
    frontCnt = 0
    for chro in population[1:]:
        #check chro is not dominated by any chro in ParetoFronts[frontCnt]
        anyDominating = False
        for frontChro in ParetoFronts[frontCnt]:
            if frontChro.fitnessesDominating(chro):
                anyDominating = True
                break
            
        if anyDominating:
            #因為point可被current ParetoFronts中某一點dominating
            #所以此點必屬於下一級的ParetoFront
            ParetoFronts.append([])
            frontCnt += 1
            ParetoFronts[frontCnt].append(chro)
            
        else:
            #find lowest front b s.t. Paretofronts[b] not dominating chro
            b = 0
            for cnt in xrange(frontCnt+1):
                if not any(frontChro2.fitnessesDominating(chro) for 
                                frontChro2 in ParetoFronts[cnt]):
                    b = cnt
                    break
            ParetoFronts[b].append(chro)
    return ParetoFronts

def highObjectivesNonDominatedSort(population):
    '''
    @param points: set for non-dominated sort
    @param M: num. of objectives or idx of objectives 
    '''
    #initialize
    for chro in population:
        chro.ParetoRank = 1
    
    #computing Pareto rank of each point
    M = population[0].numOfObjectives
    ND_helper_A(population, M)
    
    #allocating each point to corresponding front    
    numOfFronts = max(chro.ParetoRank for chro in population)
    ParetoFronts = [ [] for _ in xrange(numOfFronts)]
    [ParetoFronts[chro.ParetoRank-1].append(chro) for chro in population]
    
    return ParetoFronts

def splitSet(population):
    '''
    @return medianIdx: sorted population，第M個objective的中位數的索引值
    -如果value[medianIdx] == value[mediaIdx+1]時，medianIdx +=1，直到list結尾
    '''
    dataLen = len(population)
    M = population[0].numOfObjectives
    
    halfIdx  = (dataLen-1)/2 if dataLen % 2 else dataLen/2 -1
    #小於等於median的points分類L, 大於median的points為H
    #H當中的所有元素絕對不可能dominating L中的任一點(minimizing objective)
    for Mdx in reversed(xrange(M)):
        population.sort(key = lambda chro: chro.fitnesses[Mdx])
        values = [chro.fitnesses[Mdx] for chro in population]
        
        medianIdx = bisect.bisect_right(values, values[halfIdx])
        if medianIdx != dataLen:
            return medianIdx-1, Mdx+1
            #因為values[halfidx]一定在list內，所以medianIdx==dataLen
            #表示往後切不開，要往前切
        else:
            medianIdx = bisect.bisect_left(values, values[halfIdx])
            if medianIdx != 0:
                return medianIdx-1, Mdx+1
    #表示此集合的元素元全相同
    return -1, 0

def ND_helper_A(population, M):
    '''
    -將points依第M個objective為基準值，分為L, H兩個set
    -先求points第M個objective的中位數
    -L set中的point其第M個objective之值均小於或等於中位數
    -H set中的point其第M個objective之值均大於中位數
    '''
    dataLen = len(population)
   
    if dataLen == 2:
        #只剩下兩個點時，若能夠直接做dominating comparison，則計算值，否則不做事
        if population[0].fitnessesDominating(population[1]):
            population[1].ParetoRank = max(population[0].ParetoRank+1, population[1].ParetoRank)
        elif  population[1].fitnessesDominating(population[0]):
            population[0].ParetoRank = max(population[0].ParetoRank, population[1].ParetoRank+1)
        
    elif dataLen > 2:
        #大於兩個點時，將points依第M個objective切成L, H兩個集合
        #這邊會有問題，如果medianIdx==len(population)-1時，
        #會造成H為空集合，所以應該往下一層切
        medianIdx, newM = splitSet(population)
        
        if medianIdx == -1:
            #全部fitness值都相同，直接return
            return
        
        L, H = population[:medianIdx+1], population[medianIdx+1:]
        
        #檢查L集合那些點有dominating關係
        ND_helper_A(L, newM )       
        #在L集合中的point的Pareto rank已確定後，以此為基礎計算H中 point之Pareto rank
        ND_helper_B(L, H, newM-1)
        #檢查H集合那些點有dominating關係
        ND_helper_A(H, newM)
        

def ND_helper_B(L, H, M):
    '''
    the procedure assigns Pareto rank to the solutions in H 
    according to the solutions in L.
    假設所有L集合中的point之Pareto rank均為正確值
    
    '''
    if len(L) == 1:
        #如果L[0] dominating H中的點，則更新H中point之ParetoRank, 只有H的rank會改變
        for chro in H:
            if L[0].fitnessesDominating(chro):
                chro.ParetoRank = max(chro.ParetoRank, L[0].ParetoRank+1)
                
    elif len(H) == 1:
        #如果L中的point dominating H[0]，則更新H[0]之ParetoRank
        for chro in L:
            if chro.fitnessesDominating(H[0]):
#            if chro.dominatingM(H[0], M):
                H[0].ParetoRank = max(chro.ParetoRank+1, H[0].ParetoRank)
                
              
    elif 0<= M <= (L[0].numOfObjectives-1):
        #TODO: 在這邊改成M<=len(objectives)-1速更會比再切下去更快
        #merge, 只剩下兩個目標時，用sweep line algorithm直接求出Pareto rank
        L.sort(key=lambda chro: chro.ParetoRank, reverse=True)  #descending
        for chroH in H:
            for chroL in L:
                if chroL.fitnessesDominating(chroH):
                    chroH.ParetoRank = max(chroL.ParetoRank+1, chroH.ParetoRank)
                    break
  