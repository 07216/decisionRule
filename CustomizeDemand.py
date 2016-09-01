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
    def __init__(self,choose, t=10, limt=0, d=7, T = 20):
        self.t = t
        self.limt = min(limt,self.t)
        self.T = T
        self.d = d
        self.lenMon = 100000
        
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
    
    def avg(self,pt,pj,start,threshold,minus):
        result = 0.0
        for i in range(start,self.lenMon):
            if self.monteCarlo[pt,pj][i] > threshold:
                break
            result += self.monteCarlo[pt,pj][i]-minus
        return result,i
    
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
        
        totalLen = 1000.0
        div = int(totalLen / self.t)
        self.cons = np.zeros((self.t,2))
        for t in range(0,self.t):            
            for tt in range(t*div,(t+1)*div):
                self.cons[t,0] +=  0.25 * 1.0/totalLen * (float(tt)/totalLen) ** (6 - 1) * (1- float(tt)/totalLen) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.cons[t,1] +=  0.75 * 1.0/totalLen * (float(tt)/totalLen) ** (2 - 1) * (1- float(tt)/totalLen) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                
        self.monteCarlo = {}
        #Segmentation
        self.seg = {}
        minsup = 0.01
        mininf = 0.01
        for t in range(0,self.t):
            for j in range(0,20):
                simGamma = np.random.gamma(40,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                a = 0.
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
            
            for j in range(20,60):
                simGamma = np.random.gamma(100,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                a = 0.
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
            for j in range(0,30):
                left = 0
                leftSec = 0
                for d in range(0,self.d):
                    low = self.seg[t,2*j][d]
                    up = self.seg[t,2*j][d+1]
                    if d == 0:
                        low = 0
                    self.xi[1+(t*self.j+2*j)*self.d+d],left =  self.avg(t,2*j,left,up,low)
                    self.xi[1+(t*self.j+2*j)*self.d+d] += (up - low)* (self.lenMon - left)
                    self.xi[1+(t*self.j+2*j)*self.d+d] /= self.lenMon
                    
                    low = self.seg[t,2*j+1][d]
                    up = self.seg[t,2*j+1][d+1]
                    if d == 0:
                        low = 0
                    self.xi[1+(t*self.j+2*j+1)*self.d+d],leftSec =  self.avg(t,2*j+1,leftSec,up,low)
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] += (up - low)* (self.lenMon - leftSec)
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] /= self.lenMon
                    
        #print self.xi
        #print self.h
    
    def produceDemandForFirstCaseInResolve(self):        
        self.realDemand = []
        for t in range(0,self.t):
            for j in range(0,10):
                g = np.random.gamma(40)
                self.realDemand += [np.random.poisson(g * self.cons[t][0])]
                self.realDemand += [np.random.poisson(g * self.cons[t][1])]
            for j in range(10,30):
                g = np.random.gamma(100)
                self.realDemand += [np.random.poisson(g * self.cons[t][0])]
                self.realDemand += [np.random.poisson(g * self.cons[t][1])]
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
        
        totalLen = 1000.0
        div = int(totalLen / self.t)
        self.cons = np.zeros((self.t,2))
        for t in range(0,self.t):            
            for tt in range(t*div,(t+1)*div):
                self.cons[t,0] +=  0.25 * 1.0/totalLen * (float(tt)/totalLen) ** (6 - 1) * (1- float(tt)/totalLen) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.cons[t,1] +=  0.75 * 1.0/totalLen * (float(tt)/totalLen) ** (2 - 1) * (1- float(tt)/totalLen) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                
        self.monteCarlo = {}
        #Segmentation 
        self.h = []
        self.seg = {}
        minsup = 0.01
        mininf = 0.01
        for t in range(0,self.t):
            for j in range(0,20):
                simGamma = np.random.gamma(60,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
                
            for j in range(20,44):
                simGamma = np.random.gamma(150,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
                
            for j in range(44,60):
                simGamma = np.random.gamma(100,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                    b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                    a = self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
                
        #Expectation of arrival process
        self.xi = np.zeros((self.t*self.j*self.d+1,1),dtype=np.float)
        self.xi[0] = 1
        for t in range(0,self.t):
            for j in range(0,30):
                left = 0
                leftSec = 0
                for d in range(0,self.d):
                    low = self.seg[t,2*j][d]
                    up = self.seg[t,2*j][d+1]
                    if d == 0:
                        low = 0
                    self.xi[1+(t*self.j+2*j)*self.d+d],left =  self.avg(t,2*j,left,up,low)
                    self.xi[1+(t*self.j+2*j)*self.d+d] += (up - low)* (self.lenMon - left)
                    self.xi[1+(t*self.j+2*j)*self.d+d] /= self.lenMon
                    
                    low = self.seg[t,2*j+1][d]
                    up = self.seg[t,2*j+1][d+1]
                    if d == 0:
                        low = 0
                    self.xi[1+(t*self.j+2*j+1)*self.d+d],leftSec =  self.avg(t,2*j+1,leftSec,up,low)
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] += (up - low)* (self.lenMon - leftSec)
                    self.xi[1+(t*self.j+2*j+1)*self.d+d] /= self.lenMon 
        #Up Bound And Low Bound For Demand
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
                g = np.random.gamma(60)
                self.realDemand += [np.random.poisson(g * self.cons[t][0])]
                self.realDemand += [np.random.poisson(g * self.cons[t][1])]
            for j in range(10,22):
                g = np.random.gamma(150)
                self.realDemand += [np.random.poisson(g * self.cons[t][0])]
                self.realDemand += [np.random.poisson(g * self.cons[t][1])]
            for j in range(22,30):
                g = np.random.gamma(100)
                self.realDemand += [np.random.poisson(g * self.cons[t][0])]
                self.realDemand += [np.random.poisson(g * self.cons[t][1])]
        return self.realDemand

    def reductionALPReadIn(self):
        reader = Input.Input()
        reader.readIn('data/rm_200_4_1.6_4.0.txt')
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
        self.prob = np.zeros((self.t*self.T*self.j+1,1),dtype=np.float)
        self.c = np.zeros((self.i,1),dtype=np.float)
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
        self.prob = np.zeros((self.t*self.T,self.j), dtype=np.float)
        for t in range(0,self.t * self.T):
            for j in range(0,self.j):
                index = self.rALP.prdic[j]
                self.prob[t,j] = self.rALP.bdic[(t,index[0],index[1],index[2])]
        self.h = {}
        self.monteCarlo = {}
        self.seg = {}
        minsup = 0.001
        mininf = 0.01
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.monteCarlo[t,j] = []
                for k in range(0,self.lenMon):
                    self.monteCarlo[t,j].append(np.sum(np.random.uniform(size=self.T)<self.prob[t*self.T:(t+1)*self.T,j]))
                self.monteCarlo[t,j].sort()
                b = self.monteCarlo[t,j][int(np.ceil((1-minsup)*self.lenMon-1))]
                print b
                a = 0#self.monteCarlo[t,j][int(np.floor(mininf*self.lenMon))]
                new = []
                for d in range(0,self.d):
                    new += [float(b-a)/self.d*d+a]
                new += [b]
                self.seg[t,j] = new
                       
        #Expectation of arrival process
        self.xi = np.zeros((self.t*self.j*self.d+1,1),dtype=np.float)
        self.xi[0] = 1
        for t in range(0,self.t):
            for j in range(0,self.j):
                left = 0
                for d in range(0,self.d):
                    low = self.seg[t,j][d]
                    up = self.seg[t,j][d+1]
                    if d == 0:
                        low = 0
                    self.xi[1+(t*self.j+j)*self.d+d],left =  self.avg(t,j,left,up,low)
                    self.xi[1+(t*self.j+j)*self.d+d] += (up - low)* (self.lenMon - left)
                    self.xi[1+(t*self.j+j)*self.d+d] /= self.lenMon
        
        #construct v
        for j in range(0,self.j):
            self.v[j] = self.rALP.pval[j]
        #for t in range(0,self.t):
            #for j in range(0,self.j):
                #self.xi[1+t*self.j+j] *= 100
        #print sum(self.xi)

    def produceDemandForrALP(self):
        tmpDemand = np.zeros((self.t*self.T,self.j), dtype=np.float)
        for t in range(0,self.t*self.T):
            p = np.random.uniform()
            for j in range(0,self.j):
                p-= self.prob[t,j]
                if p<=0 :
                    tmpDemand[t,j] += 1
                    break
        self.realDemand = []
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.realDemand += [np.sum(tmpDemand[t*self.T:(t+1)*self.T,j])]
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
        