# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:42:43 2016

@author: Zhan
"""

import numpy as np
from scipy.special import gamma
from scipy.stats import gamma as Gamma
from scipy.stats import poisson

class CustomizeDemand:
    def __init__(self):
        self.sim = self.produceDemandForFirstCaseInResolve
        #self.sim = self.produceDemandForrALP
    
    def resolveDemandFirstCase(self):   
        self.i = 10
        self.j = 60
        self.t = 16
        self.limt = min(0,self.t)#observed history
        
        #Fare
        self.v = np.zeros((self.j,1))
        for j in range(0,10):
            self.v[2*j] = 300
            self.v[2*j+1] = 80
        for j in range(10,30):
            self.v[2*j] = 500
            self.v[2*j+1] = 100
            
        #Capacity
        self.c = np.zeros((self.i,1))
        for i in range(0,10):
            self.c[i] = 400
        
        #matrix A I * J
        self.A = np.zeros((self.i,self.j),dtype=np.int8)
        self.refJ = {}
        for j in range(0,10):
            self.A[j,2*j] = 1
            self.A[j,2*j+1] = 1
            self.refJ[2*j] = [j]
            self.refJ[2*j+1] = [j]
        for a in range(0,5):
            for b in range(0,5):
                if a == b:
                    continue
                j = j + 1
                self.A[a*2,2*j] = 1
                self.A[b*2+1,2*j] = 1
                self.A[a*2,2*j+1] = 1
                self.A[b*2+1,2*j+1] = 1
                self.refJ[2*j] = [a*2,b*2+1]
                self.refJ[2*j+1] = [a*2,b*2+1]
            
        #Expectation of arrival process
        self.xi = np.zeros((2*(self.t*self.j)+1,1),dtype=np.float)
        self.xi[0] = 1
        for t in range(0,self.t):
            for j in range(0,10):
                self.xi[1+t*self.j+2*j] = 40 * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.xi[1+t*self.j+2*j+1] = 40 * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(10,30):
                self.xi[1+t*self.j+2*j] = 100 * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.xi[1+t*self.j+2*j+1] = 100 * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
        
        #Up Bound And Low Bound For Demand
        self.h = np.zeros((2*(self.t*self.j+1),1))
        self.h[0] = 1
        self.h[1] = -1
        minp = 1e-2
        for t in range(0,self.t):
            for j in range(0,10):
                self.h[2*(t*self.j+1)+4*j] = Gamma.ppf(1-minp,40) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+1] = -Gamma.ppf(minp,40) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+2] = Gamma.ppf(1-minp,40) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+3] = -Gamma.ppf(minp,40) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(10,30):      
                self.h[2*(t*self.j+1)+4*j] = Gamma.ppf(1-minp,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+1] = -Gamma.ppf(minp,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+2] = Gamma.ppf(1-minp,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+3] = -Gamma.ppf(minp,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                        
        
        #print self.xi
        #print self.h
    
    def produceDemandForFirstCaseInResolve(self):        
        self.realDemand = []
        for t in range(0,self.t):
            for j in range(0,10):
                self.realDemand += [np.random.gamma(40) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)]
                self.realDemand += [np.random.gamma(40) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)]
            for j in range(10,30):
                self.realDemand += [np.random.gamma(100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)]
                self.realDemand += [np.random.gamma(100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)]
        return self.realDemand
                

    def reducationALP(self,relrALP):

        self.i = relrALP.f
        self.j = relrALP.i
        self.t = relrALP.n

        
        self.t = 20#customized periods
        self.T = 10#how a periods long
        self.limt = 2#observed history
        
        
        #Sparse P
        #self.P = np.zeros((20000,20000),dtype=np.float)
        
        #construct A
        self.A = np.zeros((self.i,self.j),dtype=np.float)
        self.xi = np.zeros((self.t*self.j+1,1),dtype=np.float)
        self.prob = np.zeros((self.t*self.T*self.j+1,1),dtype=np.float)
        self.c = np.zeros((self.i,1),dtype=np.float)
        self.h = np.zeros((2*(self.t*self.j+1),1),dtype=np.float)
        #Sparse W
        #self.W = np.zeros((2*(self.t*self.j+1),self.t*self.j+1),dtype=np.float)
        self.v = np.zeros((self.j,1),dtype=np.float)
        self.rALP = relrALP
        #look resources I required by J
        self.refJ = relrALP.dicAforJ
                #construct A
        for item in self.rALP.listA:
            #Hint Item[1] is from 1 to ...
            #print item
            self.A[item[0],item[1]-1] = 1
        #construct c
        for item in range(0,self.i):
            self.c[item] = self.rALP.flight[item][2]
            #self.c[item] = 2
        #construct Exi
        self.xi[0] = 1
        self.prob[0] = 1
        for t in range(0,self.t * self.T):
            for j in range(0,self.j):
                index = self.rALP.prdic[j]
                self.prob[1+t*self.j+j] = self.rALP.bdic[(t,index[0],index[1],index[2])]
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.xi[1+t*self.j+j] = 0
                for T in range(0,self.T):
                    self.xi[1+t*self.j+j] += self.prob[1+(t*self.T+T)*self.j+j]        
        minp = 5e-2
        #construct h
        self.h[0] = 1
        self.h[1] = -1
        for t in range(0,self.t):
            for j in range(0,self.j):                
                if self.xi[t*self.j+j+1] == 0:
                    self.h[2*(1+t*self.j+j)] = 0
                    self.h[2*(1+t*self.j+j)+1] = 0
                else:
                    self.h[2*(1+t*self.j+j)] = poisson.ppf(1-minp,self.xi[t*self.j+j+1])
                    self.h[2*(1+t*self.j+j)+1] = -poisson.ppf(minp,self.xi[t*self.j+j+1])
        #construct v
        for j in range(0,self.j):
            self.v[j] = self.rALP.pval[j]
        
        #for t in range(0,self.t):
            #for j in range(0,self.j):
                #self.xi[1+t*self.j+j] *= 100
        #print sum(self.xi)

    def produceDemandForrALP(self):
        self.realDemand = []
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.realDemand += [np.random.poisson(self.xi[1+t*self.j+j],1)]
        return self.realDemand
                
    def cal(self,t,j,minp):
        s = 0
        for T in range(0,self.T):
            index = self.rALP.prdic[j]
            p = self.rALP.bdic[(t*self.T+T,index[0],index[1],index[2])]
            i = 0
            q = p
            tmp = 1
            while(q>minp):
                i = i+1
                tmp = tmp*(i+1)
                q = (p**(i+1))/tmp
            s += i
        return s