# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 17:11:55 2016

@author: Zhan
"""

import numpy as np

class simulation:
    def __init__(self,decisionSolver,recorder):
        self.i = recorder.i
        self.j = recorder.j
        self.t = recorder.t
        self.A = recorder.A
        self.xi = recorder.xi
        self.v = recorder.v
        self.c = recorder.c
        self.h = recorder.h
        self.refJ = recorder.refJ
        self.t = recorder.t
        self.limt  = recorder.limt
        self.sim = recorder.sim
        self.d = recorder.d
        self.q = recorder.q
        self.qq = recorder.qq
        self.x = decisionSolver.x
        self.xx = decisionSolver.xx
        
    def initX(self):
        self.X = {}
        for t in range(0,self.limt):
            self.X[t] = np.zeros((self.j,1+t*self.j*self.d),dtype=np.float)
            for j in range(0,self.j):
                for p in range(0,1+t*self.j*self.d):
                    self.X[t][j,p] = self.x[t,j,p].X
                    if p!=0 and self.X[t][j,p] !=0:
                        print self.X[t][j,p],t,j,p
        for t in range(self.limt,self.t):
            self.X[t] = np.zeros((self.j,1+self.limt*self.j*self.d),dtype=np.float)
            for j in range(0,self.j):
                for p in range(0,1+self.limt*self.j*self.d):
                    self.X[t][j,p] = self.x[t,j,p].X
                    if p!=0 and self.X[t][j,p] !=0:
                        print self.X[t][j,p],t,j,p
        self.XX = {}
        for t in range(0,self.t):
            for j in range(0,self.j):
                for d in range(0,self.d):
                    self.XX[t,j,d] = self.xx[t,j,d].X
                    if self.XX[t,j,d]!=0:
                        print "XX:",self.XX[t,j,d],t,j,d
                    
    def atLeastOne(self,x):
        return int(x)+1
    
    def nonLinearDemand(self,realDemand):
        result = [0] * (self.t*self.j*self.d)
        for t in range(0,self.t):
            for j in range(0,self.j):
                for d in range(0,self.d):
                    if realDemand[t*self.j+j] <= self.h[2*(1+(t*self.j+j)*self.d+d)]:
                        result[(t*self.t+j)*self.d+d] = 1
                        break
    def aSim(self):

        realDemand = self.sim()   
        c = np.copy(self.c)
        demand = [1]
        history = np.array(demand)
        benefit = 0.0
        rplc = self.atLeastOne
        #rplc = np.ceil
        
        for t in range(0,self.limt):
            product = np.dot(self.X[t],history)
            #print product
            tmpDemand = self.nonLinearDemand(realDemand[t*self.j:(t+1)*self.j])
            for j in range(0,self.j):
                if product[j]<0:
                    print "Strange!"
                sell = max(0,min(tmpDemand[j],rplc(product[j])))
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
                            
                #benefit += int(product[j]) * self.v[j]
            demand += tmpDemand
            if self.limt != 0:
                history = np.array(demand)
        
        for t in range(self.limt,self.t):
            product = np.dot(self.X[t],history)
            #print product
            tmpDemand = self.nonLinearDemand(realDemand[t*self.j:(t+1)*self.j])
            for j in range(0,self.j):
                if product[j]<0:
                    print "Strange!"
                sell = max(0,min(tmpDemand[j],rplc(product[j])))
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
            demand = [1] + demand[(self.j+1):] + tmpDemand
            if self.limt != 0:
                history = np.array(demand)
        
        for i in range(0,self.i):
            if c[i] <0:
                print "Alert!"
                
        return benefit    

    def bookLimSim(self,x):
        realDemand = self.sim()   
        
        c = np.copy(self.c)
        benefit = 0.0
        
        for t in range(0,self.t):
            product = x[t]
            #print product
            tmpDemand = realDemand[t*self.j:(t+1)*self.j]
            for j in range(0,self.j):
                if product[j]<0:
                    print "Strange!"
                sell = max(0,min(tmpDemand[j],product[j]))
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
        
        for i in range(0,self.i):
            if c[i] <0:
                print "Alert!"
                
        return benefit    
        
    def simWithGivenDemand(self,x,realDemand):
        c = np.copy(self.c)
        benefit = 0.0
        
        for t in range(0,self.t):
            product = x[t]
            #print product
            tmpDemand = realDemand[t*self.j:(t+1)*self.j]
            for j in range(0,self.j):
                if product[j]<0:
                    print "Strange!"
                sell = max(0,min(tmpDemand[j],product[j]))
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
        
        for i in range(0,self.i):
            if c[i] <0:
                print "Alert!"
                
        return benefit    
    
    def run(self,n):
        s = 0.0
        for i in range(0,n):
            s += self.aSim()
        s /= n
        return s
    
    def bookLimRun(self,n,x):
        s = 0.0
        for i in range(0,n):
            s += self.bookLimSim(x)
        s /= n
        return s