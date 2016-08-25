# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:42:43 2016

@author: Zhan
"""

import numpy as np
from scipy.special import gamma
from scipy.stats import gamma as Gamma
from scipy.stats import poisson
from scipy.integrate import quad
import Input

class CustomizeDemand:
    def __init__(self,choose):
        self.t = 10
        self.limt = min(5,self.t)
        self.T = 100
        self.d = 7
        self.q = 3
        self.qq = self.d
        
        if choose ==0:
            self.reader = self.reductionALPReadIn()
            self.reductionALP(self.reader)
            self.sim = self.produceDemandForrALP
        elif  choose == 1:
            self.resolveDemandFirstCase()
            self.sim = self.produceDemandForFirstCaseInResolve
        elif choose == 2:
            self.resolveDemandSecondCase()
            self.sim = self.produceDemandForSecondCaseInResolve
    
    def resolveDemandFirstCase(self):   
        self.i = 10
        self.j = 60
        
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
        #Whether we need to set h sepeart from seg? h for range ,seg for segmentation
        self.h = []
        
        #Segmentation 
        self.seg = {}
        minsup = 0.01
        mininf = 0.01
        for t in range(0,self.t):
            for j in range(0,20):
                if j%2 ==0:
                    b = Gamma.ppf(1-minsup,40) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                    a = Gamma.ppf(mininf,40) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                else:
                    b = Gamma.ppf(1-minsup,40) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                    a = Gamma.ppf(mininf,40) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
                
            for j in range(20,60):
                if j%2 ==0:
                    b = Gamma.ppf(1-minsup,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                    a = Gamma.ppf(mininf,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                else:
                    b = Gamma.ppf(1-minsup,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                    a = Gamma.ppf(mininf,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
        #print self.seg
        #Expectation of arrival process
        self.xi = np.zeros((self.t*self.j*self.d+1,1),dtype=np.float)
        self.xi[0] = 1
        for t in range(0,self.t):
            s = Gamma.ppf(mininf,40)
            e = Gamma.ppf(1-minsup,40)
            k = float(e-s)/self.d
            for j in range(0,10):
                for d in range(0,self.d):
                    ss = k*d+s
                    ee = k*(d+1)+s
                    low = self.seg[t,2*j][d]
                    up = self.seg[t,2*j][d+1]
                    if d == 0:
                        low = 0
                        ss = 0
                    if d == self.d:
                        ee = Gamma.ppf(0.001,40)
                    ct = 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                    self.xi[1+(t*self.j+2*j)*self.d+d] =  quad(lambda x:(x*ct-low)*Gamma.pdf(x,40),ss,ee)[0] 
                    self.xi[1+(t*self.j+2*j)*self.d+d] += (up - low)* (1 - Gamma.cdf(ee,40))
                    
                    ss = k*d+s
                    ee = k*(d+1)+s
                    low = self.seg[t,2*j+1][d]
                    up = self.seg[t,2*j+1][d+1]
                    if d == 0:
                        low = 0
                        ss = 0
                    if d == self.d:
                        ee = Gamma.ppf(0.001,40)
                    ct = 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)                    
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] =  quad(lambda x:(x*ct-low)*Gamma.pdf(x,40),ss,ee)[0] 
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] += (up - low)* (1 - Gamma.cdf(ee,40)) 
                    
            s = Gamma.ppf(mininf,100)
            e = Gamma.ppf(1-minsup,100)
            k = float(e-s)/self.d
            for j in range(10,30):
                for d in range(0,self.d):
                    ss = k*d+s
                    ee = k*(d+1)+s
                    low = self.seg[t,2*j][d]
                    up = self.seg[t,2*j][d+1]
                    if d == 0:
                        low = 0
                        ss = 0
                    if d == self.d:
                        ee = Gamma.ppf(0.001,100)
                    ct = 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                    self.xi[1+(t*self.j+2*j)*self.d+d] =  quad(lambda x:(x*ct-low)*Gamma.pdf(x,100),ss,ee)[0] 
                    self.xi[1+(t*self.j+2*j)*self.d+d] += (up - low)* (1 - Gamma.cdf(ee,100))
                    
                    ss = k*d+s
                    ee = k*(d+1)+s
                    low = self.seg[t,2*j+1][d]
                    up = self.seg[t,2*j+1][d+1]
                    if d == 0:
                        low = 0
                        ss = 0
                    if d == self.d:
                        ee = Gamma.ppf(0.001,100)
                    ct = 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)                    
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] =  quad(lambda x:(x*ct-low)*Gamma.pdf(x,100),ss,ee)[0] 
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] += (up - low)* (1 - Gamma.cdf(ee,100)) 
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
                
    def resolveDemandSecondCase(self):
        self.i = 10
        self.j = 60
        
        #Fare
        self.v = np.zeros((self.j,1))
        for j in range(0,10):
            self.v[2*j] = 300
            self.v[2*j+1] = 80
        for j in range(10,22):
            self.v[2*j] = 500
            self.v[2*j+1] = 100
        for j in range(22,30):
            self.v[2*j] = 700
            self.v[2*j+1] = 200
            
        #Capacity
        self.c = np.zeros((self.i,1))
        for i in range(0,10):
            self.c[i] = 400
        self.c[4] = self.c[5] = 1000
        
        #matrix A I * J
        self.A = np.zeros((self.i,self.j),dtype=np.int8)
        self.refJ = {}
        for j in range(0,10):
            self.A[j,2*j] = 1
            self.A[j,2*j+1] = 1
            self.refJ[2*j] = [j]
            self.refJ[2*j+1] = [j]
        self.pushTwoLeg(0,1,0,3,1,2,0)
        self.pushTwoLeg(0,5,0,4,5,1,1)
        self.pushTwoLeg(1,5,2,4,5,3,2)
        self.pushTwoLeg(2,4,6,5,4,7,3)
        self.pushTwoLeg(2,3,6,8,9,7,4)
        self.pushTwoLeg(3,4,9,5,4,8,5)
        self.pushThreeLeg(0,2,0,4,7,6,5,1,6)
        self.pushThreeLeg(0,3,0,4,8,9,5,1,7)
        self.pushThreeLeg(1,2,2,4,7,6,5,3,8)
        self.pushThreeLeg(1,3,2,4,8,9,5,3,9)
        
        #Expectation of arrival process
        self.xi = np.zeros((2*(self.t*self.j)+1,1),dtype=np.float)
        self.xi[0] = 1
        for t in range(0,self.t):
            for j in range(0,10):
                self.xi[1+t*self.j+2*j] = 60 * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.xi[1+t*self.j+2*j+1] = 60 * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(10,22):
                self.xi[1+t*self.j+2*j] = 150 * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.xi[1+t*self.j+2*j+1] = 150 * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(22,30):
                self.xi[1+t*self.j+2*j] = 100 * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.xi[1+t*self.j+2*j+1] = 100 * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
        
        #Up Bound And Low Bound For Demand
        self.h = np.zeros((2*(self.t*self.j+1),1))
        self.h[0] = 1
        self.h[1] = -1
        minsup = 1e-1
        mininf = 7e-1
        for t in range(0,self.t):
            for j in range(0,10):
                self.h[2*(t*self.j+1)+4*j] = Gamma.ppf(1-minsup,60) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+1] = -Gamma.ppf(mininf,60) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+2] = Gamma.ppf(1-minsup,60) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+3] = -Gamma.ppf(mininf,60) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(10,22):      
                self.h[2*(t*self.j+1)+4*j] = Gamma.ppf(1-minsup,150) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+1] = -Gamma.ppf(mininf,150) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+2] = Gamma.ppf(1-minsup,150) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+3] = -Gamma.ppf(mininf,150) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
            for j in range(22,30):      
                self.h[2*(t*self.j+1)+4*j] = Gamma.ppf(1-minsup,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+1] = -Gamma.ppf(mininf,100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+2] = Gamma.ppf(1-minsup,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                self.h[2*(t*self.j+1)+4*j+3] = -Gamma.ppf(mininf,100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                        
        #print self.xi
        #print self.h
    
    def pushTwoLeg(self,s,t,a,b,c,d,id):
        self.A[a,20+4*id] = self.A[b,20+4*id] = 1
        self.A[a,20+4*id+1] = self.A[b,20+4*id+1] = 1
        self.A[c,20+4*id+2] = self.A[d,20+4*id+2] = 1
        self.A[c,20+4*id+3] = self.A[d,20+4*id+3] = 1
        self.refJ[20+4*id] = [a,b]
        self.refJ[20+4*id+1] = [a,b]
        self.refJ[20+4*id+2] = [c,d]
        self.refJ[20+4*id+3] = [c,d]
    
    def pushThreeLeg(self,s,t,a,b,c,d,e,f,id):
        self.A[a,20+4*id] = self.A[b,20+4*id] = self.A[c,20+4*id] = 1
        self.A[a,20+4*id+1] = self.A[b,20+4*id+1] = self.A[c,20+4*id+1] = 1
        self.A[d,20+4*id+2] = self.A[e,20+4*id+2] = self.A[f,20+4*id+2] = 1
        self.A[d,20+4*id+3] = self.A[e,20+4*id+3] = self.A[f,20+4*id+3] = 1
        self.refJ[20+4*id] = [a,b,c]
        self.refJ[20+4*id+1] = [a,b,c]
        self.refJ[20+4*id+2] = [d,e,f]
        self.refJ[20+4*id+3] = [d,e,f]

    def produceDemandForSecondCaseInResolve(self):        
        self.realDemand = []
        for t in range(0,self.t):
            for j in range(0,10):
                self.realDemand += [np.random.gamma(60) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)]
                self.realDemand += [np.random.gamma(60) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)]
            for j in range(10,22):
                self.realDemand += [np.random.gamma(150) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)]
                self.realDemand += [np.random.gamma(150) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)]
            for j in range(22,30):
                self.realDemand += [np.random.gamma(100) * 0.25 * 1.0/self.t * (float(t)/self.t) ** (6 - 1) * (1- float(t)/self.t) ** (2-1) * gamma(8)/gamma(2)/gamma(6)]
                self.realDemand += [np.random.gamma(100) * 0.75 * 1.0/self.t * (float(t)/self.t) ** (2 - 1) * (1- float(t)/self.t) ** (6-1) * gamma(8)/gamma(2)/gamma(6)]
        return self.realDemand

    def reductionALPReadIn(self):
        reader = Input.Input()
        reader.readIn('data/rm_200_4_1.0_4.0.txt')
        reader.product()
        reader.leg()
        reader.construct()
        return reader
        
    def reductionALP(self,relrALP):

        self.i = relrALP.f
        self.j = relrALP.i
        
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
        minsup = 1e-1
        mininf = 7e-1
        #construct h
        self.h[0] = 1
        self.h[1] = -1
        for t in range(0,self.t):
            for j in range(0,self.j):                
                if self.xi[t*self.j+j+1] == 0:
                    self.h[2*(1+t*self.j+j)] = 0
                    self.h[2*(1+t*self.j+j)+1] = 0
                else:
                    self.h[2*(1+t*self.j+j)] = poisson.ppf(1-minsup,self.xi[t*self.j+j+1])
                    self.h[2*(1+t*self.j+j)+1] = -poisson.ppf(mininf,self.xi[t*self.j+j+1])
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