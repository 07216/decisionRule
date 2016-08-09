# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 17:11:55 2016

@author: Zhan
"""

import numpy as np

class simulation:
    def __init__(self,decisionSolver):
        self.xi = decisionSolver.xi
        self.t = decisionSolver.t
        self.T = decisionSolver.T
        self.limt = decisionSolver.limt
        self.i = decisionSolver.i
        self.j = decisionSolver.j
        self.x = decisionSolver.x
        self.A = decisionSolver.A
        self.c = decisionSolver.c
        self.v = decisionSolver.v
        self.refJ = decisionSolver.refJ
    
    def construct(self):
        self.X = {}
        for t in range(0,self.limt):
            self.X[t] = np.zeros((self.j,1+t*self.j),dtype=np.float)
            for j in range(0,self.j):
                for p in range(0,1+t*self.j):
                    self.X[t][j,p] = self.x[t,j,p].X
                    if self.X[t][j,p] !=0:
                        print self.X[t][j,p],t,j,p
                    #print self.X[t]
        for t in range(self.limt,self.t):
            self.X[t] = np.zeros((self.j,1+self.limt*self.j),dtype=np.float)
            for j in range(0,self.j):
                for p in range(0,1+self.limt*self.j):
                    self.X[t][j,p] = self.x[t,j,p].X
                    if self.X[t][j,p] !=0:
                        print self.X[t][j,p],t,j,p
                    
    def calDemand(self,t):
        demand = {}
        for j in range(0,self.j):
            p = self.xi[1+t*self.j+j]
            '''
            if np.random.uniform()<p:
                demand[j] = 1
            else:
                demand[j] = 0
            '''
            demand[j] = np.random.poisson(p)
        #print demand
        return demand.values()
    
    def aSim(self):
        c = np.copy(self.c)
        demand = [1]
        history = np.array(demand)
        benefit = 0.0
        
        for t in range(0,self.limt):
            #print history
            #print self.X[t*self.period]
            product = np.dot(self.X[t],history)
            #print product
            tmpDemand = self.calDemand(t)
            for j in range(0,self.j):
                #print product[j]
                benefit += min(tmpDemand[j],int(product[j])+1) * self.v[j]
                #benefit += min(tmpDemand[j],product[j]) * self.v[j]
                #benefit += product[j] * self.v[j] *self.xi[1+t*self.j+j]
                #benefit += min(tmpDemand[j],1) * self.v[j]
                if min(tmpDemand[j],int(product[j])+1) != 0:
                    #print min(tmpDemand[j],int(product[j]))
                    for k in self.refJ[j]:
                        c[k] -= min(tmpDemand[j],int(product[j])+1)
            demand += tmpDemand
            history = np.array(demand)
            
        for t in range(self.limt,self.t):
            product = np.dot(self.X[t],history)
            #print product
            tmpDemand = self.calDemand(t)
            for j in range(0,self.j):
                #print product[j]
                benefit += min(tmpDemand[j],int(product[j])+1) * self.v[j]
                #benefit += min(tmpDemand[j],product[j]) * self.v[j]
                #benefit += product[j] * self.v[j] *self.xi[1+t*self.j+j]
                #benefit += min(tmpDemand[j],1) * self.v[j]
                if min(tmpDemand[j],int(product[j])+1) != 0:
                    #print min(tmpDemand[j],int(product[j]))
                    for k in self.refJ[j]:
                        c[k] -= min(tmpDemand[j],int(product[j])+1)
            demand = [1] + demand[(self.j+1):] + tmpDemand
            history = np.array(demand)
        #print c
        return benefit
    
    def run(self,n):
        s = 0.0
        for i in range(0,n):
            s += self.aSim()
        s /= n
        print s