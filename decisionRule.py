# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:23:11 2016

@author: Zhan
"""

from gurobipy import *
from parmap import starmap
import numpy as np

class decisionRule:
    def __init__(self):
        self.m = Model("Decision Rule Approch")
        self.kernel = 32
        
        #Sparse P
        #self.P = np.zeros((20000,20000),dtype=np.float)
        #self.prob = np.zeros((self.t*self.T*self.j+1,1),dtype=np.float)
        #Sparse W
        #self.W = np.zeros((2*(self.t*self.j+1),self.t*self.j+1),dtype=np.float)

    def inputDemand(self,recorder):
        self.i = recorder.i
        self.j = recorder.j
        self.t = recorder.t
        self.A = recorder.A
        self.xi = recorder.xi
        self.v = recorder.v
        self.c = recorder.c
        self.h = recorder.h
        self.seg = recorder.seg
        self.refJ = recorder.refJ
        self.t = recorder.t
        self.limt  = recorder.limt
        self.d = recorder.d
        self.q = recorder.q
        self.qq = recorder.qq
           
    def echoInput(self):
        print "i:%d"%self.i
        print "j:%d"%self.j
        print "t:%d"%self.t
        print "legs:"
        print self.c
        print "value:"
        print self.v
        print "h:"
        print self.h
        print "A:"
        print self.A
        print "xi:"
        print self.xi
    
    def echoVal(self):
        for x in self.h:
            print x
    
    def addVar(self):
        #add X
        self.x = {}
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j*self.d+1):
                    self.x[t,j,p] = self.m.addVar(lb=-GRB.INFINITY, name = 'X %d %d %d' % (t,j,p))
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                for p in range(0,self.limt*self.j*self.d+1):
                    self.x[t,j,p] = self.m.addVar(lb=-GRB.INFINITY, name = 'X %d %d %d' % (t,j,p))
        #add xx
        self.xx = {}
        for t in range(0,self.t):
            for j in range(0,self.j):
                for d in range(0,self.d):
                    self.xx[t,j,d] = self.m.addVar(lb=-GRB.INFINITY, name = 'XX %d %d %d' % (t,j,d))
        #add Lambda
        self.l = {}
        for p in range(0,self.i):
            for i in range(0,self.t*self.j*(self.d+1)+2):
                self.l[p,i] = self.m.addVar(lb=-GRB.INFINITY, ub=0, name = 'Lambda %d %d' % (p,i))
                
        #add Gamma
        self.g = {}
        for t in range(0,self.t):
            for p in range(0,self.j):
                for i in range(0,self.t*self.j*(self.d+1)+2):
                    self.g[t,p,i] = self.m.addVar(name = 'Omega %d %d %d' %(t,p,i))
        #add Omega
        self.o = {}
        for t in range(0,self.t):
            for p in range(0,self.j):
                for i in range(0,self.t*self.j*(self.d+1)+2):
                    self.o[t,p,i] = self.m.addVar(lb=-GRB.INFINITY, ub=0, name = 'Omega %d %d %d' %(t,p,i))

        self.m.update()
        
        #x size for first limt t, (self.j , t*self.j+1)
    
    def addOpt(self):
        obj  = LinExpr()
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j*self.d+1):
                    if self.v[j,0] * self.xi[p,0] ==0:
                        continue
                    obj += self.v[j,0] * self.x[t,j,p] * self.xi[p,0]
                    #obj += self.v[j,0] * self.x[j,p]
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                #first element of xi = 1
                obj += self.v[j,0] * self.x[t,j,0]
                for p in range(1,self.limt*self.j*self.d+1):
                    if self.v[j,0] * self.xi[(t-self.limt)*self.j*self.d+p,0] ==0:
                        continue
                    obj += self.v[j,0] * self.x[t,j,p] * self.xi[(t-self.limt)*self.j*self.d+p,0]
                    
        for t in range(0,self.t):
            for j in range(0,self.j):
                for d in range(0,self.d):
                    obj += self.v[j,0] * self.xx[t,j,d] * self.xi[1+(t*self.j+j)*self.d+d,0]
                    
        self.m.setObjective(obj,GRB.MAXIMIZE)
        
    def paraLambda(self,(pt,pj,pd)):
        base = 2+(pt*self.j+pj)*(self.d+1)+pd
        lhs = {}
        rhs = {}
        for i in range(0,self.i):
            lhs[i] = LinExpr()
            rhs[i] = LinExpr()
        if self.seg[pt,pj][1] != self.seg[pt,pj][0]:
            for i in range(0,self.i):
                lhs[i] += self.l[i,base] * -1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
                lhs[i] += self.l[i,base+1] * 1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
        else:
            for i in range(0,self.i):
                lhs[i] += self.l[i,base] * -1.0 
                lhs[i] += self.l[i,base+1] * 1.0                        
        #Ax xi
        col = 1+(pt*self.j+pj)*self.d+pd
        for t in range(pt+1,self.limt):
            for j in range(0,self.j):
                for i in self.refJ[j]:
                    rhs[i] += self.x[t,j,col]
        for t in range(max(pt+1,self.limt),min(self.t,pt+1+self.limt)):
            for j in range(0,self.j):
                for i in self.refJ[j]:
                    rhs[i] += self.x[t,j,1+((pt-(t-self.limt))*self.j+pj)*self.d+pd]
        #Axx xi
        for i in self.refJ[pj]:
            rhs[i] += self.xx[pt,pj,pd]
        for i in range(0,self.i):
            self.m.addConstr(rhs[i],GRB.EQUAL,lhs[i],'Z1 %d %d %d %d' %(pt,pj,pd,i))
        return 1
        
    def paraGamma(self,(t,pt,pj,pd)):
        base = 2+(pt*self.j+pj)*(self.d+1)+pd
        lhs = {}
        rhs = {}
        for i in range(0,self.j):
            lhs[i] = LinExpr()
            rhs[i] = LinExpr()
        if self.seg[pt,pj][1] != self.seg[pt,pj][0]:
            for i in range(0,self.j):
                lhs[i] += self.g[t,i,base] * -1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
                lhs[i] += self.g[t,i,base+1] * 1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
        else:
            for i in range(0,self.j):
                lhs[i] += self.g[t,i,base] * -1.0 
                lhs[i] += self.g[t,i,base+1] * 1.0                        
        #Ax xi
        col = 1+(pt*self.j+pj)*self.d+pd
        if t < self.limt:
            if pt < t:                              
                for j in range(0,self.j):
                    rhs[j] += self.x[t,j,col]
        else:
            if pt < t and pt >= t - self.limt:                            
                for j in range(0,self.j):
                    rhs[j] += self.x[t,j,1+((pt-(t-self.limt))*self.j+pj)*self.d+pd]                            
        #Axx xi
        if pt == t:
            rhs[pj] += self.xx[pt,pj,pd]
        for j in range(0,self.j):
            self.m.addConstr(rhs[j],GRB.EQUAL,lhs[j])
        return 1
            
    def paraOmega(self,(t,pt,pj,pd)):
        base = 2+(pt*self.j+pj)*(self.d+1)+pd
        lhs = {}
        rhs = {}
        for i in range(0,self.j):
            lhs[i] = LinExpr()
            rhs[i] = LinExpr()
        if self.seg[pt,pj][1] != self.seg[pt,pj][0]:
            for i in range(0,self.j):
                lhs[i] += self.o[t,i,base] * -1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
                lhs[i] += self.o[t,i,base+1] * 1.0 / (self.seg[pt,pj][pd+1] - self.seg[pt,pj][pd])
        else:
            for i in range(0,self.j):
                lhs[i] += self.o[t,i,base] * -1.0 
                lhs[i] += self.o[t,i,base+1] * 1.0                        
        #Ax xi
        col = 1+(pt*self.j+pj)*self.d+pd
        if t < self.limt:
            if pt < t:                              
                for j in range(0,self.j):
                    rhs[j] += self.x[t,j,col]
        else:
            if pt < t and pt >= t - self.limt:                            
                for j in range(0,self.j):
                    rhs[j] += self.x[t,j,1+((pt-(t-self.limt))*self.j+pj)*self.d+pd]                            
        #Axx xi
        if pt == t:
            rhs[pj] += self.xx[pt,pj,pd] - 1
        for j in range(0,self.j):
            self.m.addConstr(rhs[j],GRB.EQUAL,lhs[j])
        return 1
            
    def addConstr(self):
        #Lambda h <= 0 
        lhs = {}
        for i in range(0,self.i):
            lhs[i] = LinExpr()
            lhs[i] += self.l[i,0] - self.l[i,1]
        #For Constant demand xi
        for pt in range(0,self.t):
            for pj in range(0,self.j):
                if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                    base = (pt*self.j+pj)*(self.d+1)+2
                    for i in range(0,self.i):
                        lhs[i] += self.l[i,base] * - float(self.seg[pt,pj][0])/self.d
                        lhs[i] += self.l[i,base+self.d] * float(self.seg[pt,pj][0])/self.d
        for i in range(0,self.i):
            self.m.addConstr(lhs[i],GRB.LESS_EQUAL,0)
        #Lambda for Lambda*W = Z1
        #Frist Column
        lhs = {}
        rhs = {}
        for i in range(0,self.i):
            lhs[i] = LinExpr()
            rhs[i] = LinExpr()
        for i in range(0,self.i):
            lhs[i] += self.l[i,0] - self.l[i,1]
        for pt in range(0,self.t):
            for pj in range(0,self.j):
                if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                    continue
                base = (pt*self.j+pj)*(self.d+1)+2
                for i in range(0,self.i):
                    lhs[i] += self.l[i,base] * float(self.seg[pt,pj][1]) / (self.seg[pt,pj][1] - self.seg[pt,pj][0])
                    lhs[i] += self.l[i,base+1] * float(self.seg[pt,pj][0]) / (-self.seg[pt,pj][1] + self.seg[pt,pj][0])
        for t in range(0,self.t):
            for j in range(0,self.j):
                for i in self.refJ[j]:
                    rhs[i] += self.x[t,j,0]
        for i in range(0,self.i):
            rhs[i] += -self.c[i,0]
        for i in range(0,self.i):
            self.m.addConstr(lhs[i], GRB.EQUAL, rhs[i],'Constant %d' %(i))
        #Beside first column
        Parallel(n_jobs=self.kernel)(delayed(self.paraLambda)(pt,pj,pd) for pt in range(self.t) for pj in range(self.j) for pd in range(self.d))
        #Gamma h >=0
        for t in range(0,self.t):
            lhs = {}
            for i in range(0,self.j):
                lhs[i] = LinExpr()
                lhs[i] += self.g[t,i,0] - self.g[t,i,1]
            #For Constant demand xi
            for pt in range(0,self.t):
                for pj in range(0,self.j):
                    if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                        base = (pt*self.j+pj)*(self.d+1)+2
                        for i in range(0,self.j):
                            lhs[i] += self.g[t,i,base] * - float(self.seg[pt,pj][0])/self.d
                            lhs[i] += self.g[t,i,base+self.d] * float(self.seg[pt,pj][0])/self.d
            for i in range(0,self.j):
                self.m.addConstr(lhs[i],GRB.GREATER_EQUAL,0)
        for t in range(0,self.t):
            lhs = {}
            rhs = {}
            for i in range(0,self.j):
                lhs[i] = LinExpr()
                rhs[i] = LinExpr()
            for j in range(0,self.j):
                lhs[j] += self.g[t,j,0] - self.g[t,j,1]
            for pt in range(0,self.t):
                for pj in range(0,self.j):
                    if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                        continue
                    base = (pt*self.j+pj)*(self.d+1)+2
                    for i in range(0,self.j):
                        lhs[i] += self.g[t,i,base] * float(self.seg[pt,pj][1]) / (self.seg[pt,pj][1] - self.seg[pt,pj][0])
                        lhs[i] += self.g[t,i,base+1] * float(self.seg[pt,pj][0]) / (-self.seg[pt,pj][1] + self.seg[pt,pj][0])
            for j in range(0,self.j):
                rhs[j] += self.x[t,j,0]
            for i in range(0,self.j):
                self.m.addConstr(lhs[i], GRB.EQUAL, rhs[i])
            #Beside first column
            Parallel(n_jobs=self.kernel)(delayed(self.paraGamma)(t,pt,pj,pd) for pt in range(self.t) for pj in range(self.j) for pd in range(self.d))
        #Omega h <=0
        for t in range(0,self.t):
            lhs = {}
            for i in range(0,self.j):
                lhs[i] = LinExpr()
                lhs[i] += self.o[t,i,0] - self.o[t,i,1]
            #For Constant demand xi
            for pt in range(0,self.t):
                for pj in range(0,self.j):
                    if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                        base = (pt*self.j+pj)*(self.d+1)+2
                        for i in range(0,self.j):
                            lhs[i] += self.o[t,i,base] * - float(self.seg[pt,pj][0])/self.d
                            lhs[i] += self.o[t,i,base+self.d] * float(self.seg[pt,pj][0])/self.d
            for i in range(0,self.j):
                self.m.addConstr(lhs[i],GRB.LESS_EQUAL,0)
        #Omega for Omega W= Z2
        #before limt
        for t in range(0,self.t):
            lhs = {}
            rhs = {}
            for i in range(0,self.j):
                lhs[i] = LinExpr()
                rhs[i] = LinExpr()
            for j in range(0,self.j):
                lhs[j] += self.o[t,j,0] - self.o[t,j,1]
            for pt in range(0,self.t):
                for pj in range(0,self.j):
                    if self.seg[pt,pj][1] == self.seg[pt,pj][0]:
                        continue
                    base = (pt*self.j+pj)*(self.d+1)+2
                    for i in range(0,self.j):
                        lhs[i] += self.o[t,i,base] * float(self.seg[pt,pj][1]) / (self.seg[pt,pj][1] - self.seg[pt,pj][0])
                        lhs[i] += self.o[t,i,base+1] * float(self.seg[pt,pj][0]) / (-self.seg[pt,pj][1] + self.seg[pt,pj][0])
            for j in range(0,self.j):
                rhs[j] += self.x[t,j,0]
            for i in range(0,self.j):
                self.m.addConstr(lhs[i], GRB.EQUAL, rhs[i])
            #Beside first column
            Parallel(n_jobs=self.kernel)(delayed(self.paraOmega)(t,pt,pj,pd) for pt in range(self.t) for pj in range(self.j) for pd in range(self.d))
                        
    def solve(self):
        self.m.optimize()
    
    def echoOpt(self):
        
        obj  = 0
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j+1):
                    if self.xi[1+t*self.j+j,0] * self.v[j,0] * self.xi[p,0] ==0:
                        continue
                    obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,p].X * self.xi[p,0] 
                    #obj += self.v[j,0] * self.x[j,p]
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                #first element of xi = 1
                obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,0].X
                for p in range(1,self.limt*self.j+1):
                    if self.xi[1+t*self.j+j,0] * self.v[j,0] * self.xi[(t-self.limt)*self.j+p,0] ==0:
                        continue
                    obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,p].X * self.xi[(t-self.limt)*self.j+p,0]
        print obj
    
    def writeMPS(self):
        self.m.write("out.lp")
    
#    def echoSolution(self):
