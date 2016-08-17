# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:57:01 2016

@author: Zhan
"""

import numpy as np
from pyevolve import *

class gene:
    def __init__(self,decisionSolver,simulator):
        self.i = decisionSolver.i
        self.j = decisionSolver.j
        self.t = decisionSolver.t
        self.x = decisionSolver.x
        self.eval = simulator.bookLimSim
        self.evalPro = simulator.bookLimRun
        self.evalWithGivenDemand = simulator.simWithGivenDemand
        
    def initX(self):
        self.X = {}
        for t in range(0,self.t):
            self.X[t] = []
            for j in range(0,self.j):
                    self.X[t] += [int(self.x[t,j,0].X)]
                    
    def setSTD(self,realDemand):
        self.realDemand = realDemand
    
    def zeroX(self):
        self.X = {}
        for t in range(0,self.t):
            self.X[t] = []
            for j in range(0,self.j):
                    self.X[t] += [0]

    def trans(self,x):
        X = {}
        for t in range(0,self.t):
            for j in range(0,self.j):
                X[t] = x[t*self.j:(t+1)*self.j]
        return X
    
    def evalPlus(self,x):
        s = 0.0
        i = 0
        for item in self.realDemand.values():
            s += self.evalWithGivenDemand(x,item)
            i += 1
        return s/i
        
    def evolveScore(self,x):
        xx = self.trans(x)
        X = self.X.copy()
        for t in range(0,self.t):
            for j in range(0,self.j):
                X[t][j] += xx[t][j]
                if X[t][j] <0:
                    X[t][j] = 0
        return self.evalPlus(X)#Run 10 times with constant simulated demand
        #return self.evalPro(10,X)#Run 10 times with changing simulated demand
        #return self.eval(X)#Run one time 
    
    def evolve(self):
        genome = G1DList.G1DList(self.t*self.j)
        genome.evaluator.set(self.evolveScore)
        genome.setParams(rangemin=0,rangemax=2)
        ga = GSimpleGA.GSimpleGA(genome)
        ga.evolve(freq_stats=10)        
        
        xx = self.trans(ga.bestIndividual())
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.X[t][j] += xx[t][j]
                if self.X[t][j] <0:
                    self.X[t][j] =0