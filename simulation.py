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
        self.refJ = recorder.refJ
        self.t = recorder.t
        self.limt  = recorder.limt
        self.sim = recorder.sim
        self.d = recorder.d
        self.fun = recorder.pieceWiseFunctionOnMesh
        self.mesh = recorder.mesh
        self.r = recorder.r
        self.a = recorder.a
        self.b = recorder.b
        self.w = recorder.w
        self.h = recorder.h
        self.x = decisionSolver.x
        self.test = np.zeros((self.t,self.j,self.a,self.b))
                        
    def bookLimLeft(self):  
        realDemand = self.sim()   
        realNonLinearDemand = self.nonLinearDemand(realDemand)
        c = np.copy(self.c)
        benefit = 0.0
        #rplc = self.identity
        #rplc = np.ceil
        #rplc = np.round
        product = {}
        cu = [0] * self.j
        for t in range(self.t):
            for j in range(self.j):
                product[t,j] = 0
                for a in range(self.a):
                    for b in range(self.b):
                        product[t,j] += self.x[t,j,a,b].x * realNonLinearDemand[self.t*self.j + (t*self.j+j)*self.r +\
                        a*self.b+b ]
                #product[t,j] = self.inter(t,j,cu[j]+(t==0))
                #cu[j] += realDemand[t*self.j+j]
                    
        for t in range(self.t):
            #print product
            productDemand = realDemand[t*self.j:(t+1)*self.j]
            for j in range(0,self.j):
                sell = product[t,j]
                for k in self.refJ[j]:
                    c[k] -= sell
        print c
    
    def echoXX(self):
        for t in range(0,self.t):
            for j in range(0,self.j):
                print self.XX[t,j]
        
    def atLeastOne(self,x):
        return int(x)+1
        
    def atLeastTwo(self,x):
        return int(x)+2

    def atLeastThree(self,x):
        return int(x)+3
    
    def nonLinearDemand(self,realDemand):
        result = realDemand
        cu = [0] * self.j
        for t in range(self.t):
            for j in range(self.j):
                add = [0] * self.r
                tmp = self.fun(t,j,realDemand[t*self.j+j],cu[j]+(t==0))
                for (a,b,value) in tmp:
                    add[a*self.b+b] += value
                cu[j] += realDemand[t*self.j+j]
                result += add
        return result

    def identity(self,x):
        return x
    
    def aSim(self):
        realDemand = self.sim()   
        realNonLinearDemand = self.nonLinearDemand(realDemand)
        c = np.copy(self.c)
        benefit = 0.0
        #rplc = self.identity
        #rplc = np.ceil
        #rplc = np.round
        product = {}
        cu = [0] * self.j
        for t in range(self.t):
            for j in range(self.j):
                product[t,j] = 0
                for a in range(self.a):
                    for b in range(self.b):
                        product[t,j] += self.x[t,j,a,b].x * realNonLinearDemand[self.t*self.j + (t*self.j+j)*self.r +\
                        a*self.b+b ]
                #product[t,j] = self.inter(t,j,cu[j]+(t==0))
                #cu[j] += realDemand[t*self.j+j]
                    
        for t in range(self.t):
            #print product
            productDemand = realDemand[t*self.j:(t+1)*self.j]
            for j in range(0,self.j):
                sell = product[t,j]
                #sell = max(0,product[j])
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
        
        for i in range(0,self.i):
            if c[i] <0:
                print "Alert!"
                
        for t in range(self.t):
            for j in range(self.j):
                for a in range(self.a):
                    for b in range(self.b):
                        self.test[t,j,a,b] += realNonLinearDemand[self.t*self.j + (t*self.j+j)*self.r +\
                        a*self.b+b ]
                        
        #if lessZero == 1:
         #   print "Strange"
        return benefit    

    def inter(self,t,j,y):
        y = max(self.mesh[t,j][1][0],min(self.mesh[t,j][1][-1],y))
        for i in range(self.b):
            if self.mesh[t,j][1][i] >= y:
                break
        if i==0:
            return self.x[t,j,self.a-1,0].x
            
        return self.x[t,j,self.a-1,i-1].x+(self.x[t,j,self.a-1,i].x-self.x[t,j,self.a-1,i-1].x)* \
                (y-self.mesh[t,j][1][i-1]) / float(self.mesh[t,j][1][i]-self.mesh[t,j][1][i-1])
        
    def bookLimSim(self,rplc):
        realDemand = self.sim()   
        c = np.copy(self.c)
        benefit = 0.0
        #rplc = self.identity
        #rplc = np.ceil
        #rplc = np.round
        product = {}
        cu = [0] * self.j
        for t in range(self.t):
            for j in range(self.j):
                product[t,j] = self.inter(t,j,cu[j]+(t==0))
                cu[j] += realDemand[t*self.j+j]
                    
        for t in range(self.t):
            #print product
            productDemand = realDemand[t*self.j:(t+1)*self.j]
            for j in range(0,self.j):
                sell = max(0,min(productDemand[j],rplc(product[t,j])))
                #sell = max(0,product[j])
                if sell != 0:
                    for k in self.refJ[j]:
                        sell = min(sell,c[k])
                    benefit += sell * self.v[j]
                    for k in self.refJ[j]:
                        c[k] -= sell
        
        for i in range(0,self.i):
            if c[i] <0:
                print "Alert!"
        #if lessZero == 1:
         #   print "Strange"
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
        '''
        for t in range(self.t):
            for j in range(self.j):
                for a in range(self.a):
                    for b in range(self.b):
                        print self.test[t,j,a,b]/float(n)-self.xi[t,j,a,b]
        '''
        return s
    
    def bookLimRun(self,rplc,n):
        s = 0.0
        for i in range(0,n):
            s += self.bookLimSim(rplc)
        s /= n
        return s