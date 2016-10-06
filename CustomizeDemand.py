# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:42:43 2016

@author: Zhan
"""

import numpy as np
import scipy as sp
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
        self.lenMon = 10000
        self.lenMonSec = 10000
        self.a = 7
        self.b = 10
        self.r = self.a * self.b
        
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
    
    def avg(self,pt,pj,start,threshold,minus,diff):
        result = 0.0
        for i in range(start,self.lenMon):
            if self.monteCarlo[pt,pj][i] > threshold:
                return result,i
            result += self.monteCarlo[pt,pj][i]-minus
        return result,self.lenMon
    
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
        
        totalLen = 1000.0
        div = int(totalLen / self.t)
        self.cons = np.zeros((self.t,2))
        for t in range(0,self.t):            
            for tt in range(t*div,(t+1)*div):
                self.cons[t,0] +=  0.25 * 1.0/totalLen * (float(tt)/totalLen) ** (6 - 1) * (1- float(tt)/totalLen) ** (2-1) * gamma(8)/gamma(2)/gamma(6)
                self.cons[t,1] +=  0.75 * 1.0/totalLen * (float(tt)/totalLen) ** (2 - 1) * (1- float(tt)/totalLen) ** (6-1) * gamma(8)/gamma(2)/gamma(6)
                
        minf = 0.01
        msup = 0.99
        mxinf = 0.01
        mxsup = 0.99
        self.monteCarlo = {}
        self.lb={}
        self.ub={}
        self.xub = {}
        self.xlb = {}
        for t in range(0,self.t):
            for j in range(0,20):
                simGamma = np.random.gamma(40,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                self.ub[t,j] = self.monteCarlo[t,j][int(np.ceil(msup*self.lenMon-1))]
                self.lb[t,j] = self.monteCarlo[t,j][int(np.floor(minf*self.lenMon))]
                self.xub[t,j] = self.monteCarlo[t,j][int(np.ceil(mxsup*self.lenMon-1))]
                self.xlb[t,j] = self.monteCarlo[t,j][int(np.floor(mxinf*self.lenMon))]
            
            for j in range(20,60):
                simGamma = np.random.gamma(100,size=(self.lenMon))
                if j%2 ==0:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    self.monteCarlo[t,j].sort()
                else:
                    self.monteCarlo[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                    self.monteCarlo[t,j].sort()
                self.ub[t,j] = self.monteCarlo[t,j][int(np.ceil(msup*self.lenMon-1))]
                self.lb[t,j] = self.monteCarlo[t,j][int(np.floor(minf*self.lenMon))]
                self.xub[t,j] = self.monteCarlo[t,j][int(np.ceil(mxsup*self.lenMon-1))]
                self.xlb[t,j] = self.monteCarlo[t,j][int(np.floor(mxinf*self.lenMon))]
                    
        self.yub = np.zeros((self.t,self.j),dtype=np.float)
        self.ylb = np.zeros((self.t,self.j),dtype=np.float)
        for t in range(1,self.t):
            for j in range(self.j):
                self.yub[t,j] = self.ub[t-1,j] + self.yub[t-1,j]
                self.ylb[t,j] = self.lb[t-1,j] + self.ylb[t-1,j]
        for j in range(self.j):
            self.yub[0,j] = 1
        self.produceMesh()
        tmp = {}
        cu = {}
        self.xi = np.zeros((self.t,self.j,self.a,self.b),dtype=np.float)
        for s in range(self.lenMonSec):
            #simulation
            for t in range(self.t):
                
                simGamma = np.random.gamma(40)
                for j in range(20):
                    if j%2 ==0:
                        tmp[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    else:
                        tmp[t,j] = np.random.poisson(simGamma * self.cons[t,1])
                
                simGamma = np.random.gamma(100)
                for j in range(20,60):
                    if j%2 ==0:
                        tmp[t,j] = np.random.poisson(simGamma * self.cons[t,0])
                    else:
                        tmp[t,j] = np.random.poisson(simGamma * self.cons[t,1])
            #calculate
            for j in range(self.j):
                cu[j] = 0
            for t in range(self.t):
                for j in range(self.j):
                    result = self.pieceWiseFunctionOnMesh(t,j,tmp[t,j],cu[j]+(t==0))
                    for (a,b,value) in result:
                        self.xi[t,j,a,b] += value
                    cu[j] += tmp[t,j]
        #mean
        for t in range(self.t):
            for j in range(self.j):
                for a in range(self.a):
                    for b in range(self.b):
                        self.xi[t,j,a,b] /= float(self.lenMonSec)
        self.buildWh()
        
    def buildWh(self):
        #W
        #First Block
        self.row = 8*self.t*self.j+2+self.t*self.j*self.r
        self.col = self.t*self.j+self.r*self.t*self.j+1
        row = []
        col = []
        data = []
        for t in range(self.t):
            for j in range(self.j):
                row += [2*(t*self.j+j),2*(t*self.j+j)+1]
                col += [t*self.j+j,t*self.j+j]
                data += [1,-1]
        #Second Block
        rowBase = self.t*self.j*2
        colBase = self.t*self.j
        for t in range(self.t):
            for j in range(self.j):
                for r in range(self.r):
                    row += [rowBase + 2*(t*self.j+j)]
                    col += [colBase + (t*self.j+j)*self.r+r]
                    data += [1]
                for r in range(self.r):
                    row += [rowBase + 2*(t*self.j+j) + 1]
                    col += [colBase + (t*self.j+j)*self.r+r]
                    data += [-1]
        #Last Block
        rowBase = 4*self.t*self.j
        for j in range(self.j):
            row += [rowBase + 4*j,rowBase + 4*j+1,rowBase + 4 * j + 2 , rowBase + 4*j+3 ]
            col += [j,j,self.t*self.j+self.t*self.j*self.r,self.t*self.j+self.t*self.j*self.r]
            data += [1,-1,1,-1]
            for a in range(self.a):
                for b in range(self.b):
                    row += [rowBase + 4*j]
                    col += [self.t*self.j+j*self.r+a*self.b+b]
                    data += [-self.mesh[0,j][0][a]]
            
            for a in range(self.a):
                for b in range(self.b):
                    row += [rowBase + 4*j+1]
                    col += [self.t*self.j+j*self.r+a*self.b+b]
                    data += [self.mesh[0,j][0][a]]
            
            for a in range(self.a):
                for b in range(self.b):
                    row += [rowBase + 4*j+2]
                    col += [self.t*self.j+j*self.r+a*self.b+b]
                    data += [-self.mesh[0,j][1][b]]
            
            for a in range(self.a):
                for b in range(self.b):
                    row += [rowBase + 4*j+3]
                    col += [self.t*self.j+j*self.r+a*self.b+b]
                    data += [self.mesh[0,j][1][b]]
        colBase = self.t*self.j
        for t in range(1,self.t):
            for j in range(self.j):
                
                row += [rowBase + 4*(t*self.j+j),rowBase + 4*(t*self.j+j)+1]
                col += [t*self.j+j,t*self.j+j]                
                data += [1,-1]
                
                for tt in range(t):
                    row += [rowBase + 4*(t*self.j+j)+2]
                    col += [tt*self.j+j]
                    data += [1]
                
                for tt in range(t):
                    row += [rowBase + 4*(t*self.j+j)+3]
                    col += [tt*self.j+j]
                    data += [-1]
                        
                for a in range(self.a):
                    for b in range(self.b):
                        row += [rowBase + 4*(t*self.j+j)]
                        col += [colBase + (t*self.j+j)*self.r+a*self.b+b]
                        data += [-self.mesh[t,j][0][a]]
                
                for a in range(self.a):
                    for b in range(self.b):
                        row += [rowBase + 4*(t*self.j+j)+1]
                        col += [colBase + (t*self.j+j)*self.r+a*self.b+b]
                        data += [self.mesh[t,j][0][a]]
                
                for a in range(self.a):
                    for b in range(self.b):
                        row += [rowBase + 4*(t*self.j+j)+2]
                        col += [colBase + (t*self.j+j)*self.r+a*self.b+b]
                        data += [-self.mesh[t,j][1][b]]
                
                for a in range(self.a):
                    for b in range(self.b):
                        row += [rowBase + 4*(t*self.j+j)+3]
                        col += [colBase + (t*self.j+j)*self.r+a*self.b+b]
                        data += [self.mesh[t,j][1][b]]
            
        row += [self.t*self.j*8,self.t*self.j*8+1]
        col += [self.t*self.j+self.t*self.j*self.r,self.t*self.j+self.t*self.j*self.r]
        data += [1,-1]
        rowBase = self.t*self.j*8+2
        colBase = self.t*self.j
        for t in range(self.t):
            for j in range(self.j):
                for  r in range(self.r):
                    row += [rowBase + (t*self.j+j)*self.r+r]
                    col += [colBase + (t*self.j+j)*self.r+r]
                    data += [1]
        
        self.w = sp.sparse.csc_matrix((data,(row,col)),\
        shape=(self.t*self.j*8+2+self.t*self.j*self.r,self.t*self.j+self.t*self.j*self.r+1))
        
        data = []
        for t in range(self.t):
            for j in range(self.j):
                data += [self.xlb[t,j]]
                data += [-self.xub[t,j]]
        for t in range(self.t):
            for j in range(self.j):
                data += [1,-1]
        for t in range(self.t):
            for j in range(self.j):
                data += [0,0,0,0]
        data += [1,-1]
        data += [0] * (self.t*self.j*self.r)
        self.h = np.array(data)

        
    def tri(self,t,j,x1,y1,x2,y2,x3,y3,x,y):#
        x1 = self.mesh[t,j][0][x1]
        x2 = self.mesh[t,j][0][x2]
        x3 = self.mesh[t,j][0][x3]
        y1 = self.mesh[t,j][1][y1]
        y2 = self.mesh[t,j][1][y2]
        y3 = self.mesh[t,j][1][y3]
        dx1 = x2-x1
        dy1 = y2-y1
        dx2 = x3-x1
        dy2 = y3-y1
        A = np.array([[dx1,dy1],[dx2,dy2]])
        b = np.array([1,1])
        xx = np.linalg.solve(A,b)
        a = xx[0]
        b = xx[1]
        d = -1
        return -(a*(x-x1)+b*(y-y1)+d)
    
    def pieceWiseFunctionOnMesh(self,t,j,x,y):#
        x = max(min(self.ub[t,j],x),self.lb[t,j])
        y = max(min(self.yub[t,j],y),self.ylb[t,j])
        flag = 0
        for a in range(self.a):
            if self.mesh[t,j][0][a] >= x:
                if self.mesh[t,j][0][a] == x:
                    flag = 1
                break
        for b in range(self.b):
            if self.mesh[t,j][1][b] >= y :
                if self.mesh[t,j][1][b] == y:
                    if flag == 1:
                        flag = 3
                    else:
                        flag = 2
                break
        if flag == 0 :
            dy = y-self.mesh[t,j][1][b-1]
            dx = x-self.mesh[t,j][0][a]
            ddy = self.mesh[t,j][1][b] - self.mesh[t,j][1][b-1]
            ddx = self.mesh[t,j][0][a-1] - self.mesh[t,j][0][a]
            if dx*ddy - ddx*dy>0:
                result = [(a-1,b,self.tri(t,j,a-1,b,a,b,a,b-1,x,y)), (a,b,self.tri(t,j,a,b,a-1,b,a,b-1,x,y)), \
                (a,b-1,self.tri(t,j,a,b-1,a,b,a-1,b,x,y))]
            else:
                result = [(a-1,b-1,self.tri(t,j,a-1,b-1,a-1,b,a,b-1,x,y)), (a-1,b,self.tri(t,j,a-1,b,a,b-1,a-1,b-1,x,y)), \
                (a,b-1,self.tri(t,j,a,b-1,a-1,b-1,a-1,b,x,y))]
        elif flag == 1:
                result = [(a,b-1,1-(y-self.mesh[t,j][1][b-1])/float(self.mesh[t,j][1][b]-self.mesh[t,j][1][b-1])), \
                        (a,b,(y-self.mesh[t,j][1][b-1])/float(self.mesh[t,j][1][b]-self.mesh[t,j][1][b-1]))]
        elif flag == 2:
                result = [(a-1,b,1-(x-self.mesh[t,j][0][a-1])/float(self.mesh[t,j][0][a]-self.mesh[t,j][0][a-1])), \
                        (a,b,(x-self.mesh[t,j][0][a-1])/float(self.mesh[t,j][0][a]-self.mesh[t,j][0][a-1]))]
        else:
            result = [(a,b,1)]
        return result
    
    def produceMesh(self):
        self.mesh = {}
        for t in range(self.t):
            for j in range(self.j):
                self.mesh[t,j] = [np.linspace(self.lb[t,j],self.ub[t,j],self.a),\
                np.linspace(self.ylb[t,j],self.yub[t,j],self.b)]
    
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
        self.monteCarlo = {}
        minf = 0.01
        msup = 0.99
        mxsup = 0.9
        mxinf = 0.1
        self.xub = np.zeros((self.t,self.j),dtype=np.float)
        self.xlb = np.zeros((self.t,self.j),dtype=np.float)
        self.ub = np.zeros((self.t,self.j),dtype=np.float)
        self.lb = np.zeros((self.t,self.j),dtype=np.float)
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.monteCarlo[t,j] = []
                for k in range(0,self.lenMon):
                    self.monteCarlo[t,j].append(np.sum(np.random.uniform(size=self.T)<self.prob[t*self.T:(t+1)*self.T,j]))
                self.monteCarlo[t,j].sort()
                self.xub[t,j] = self.monteCarlo[t,j][int(np.ceil(mxsup*self.lenMon-1))]
                self.xlb[t,j] = self.monteCarlo[t,j][int(np.floor(mxinf*self.lenMon))]
                self.ub[t,j] = self.monteCarlo[t,j][int(np.ceil(msup*self.lenMon-1))]
                #print b
                self.lb[t,j] = self.monteCarlo[t,j][int(np.floor(minf*self.lenMon))]
                       
        self.yub = np.zeros((self.t,self.j),dtype=np.float)
        self.ylb = np.zeros((self.t,self.j),dtype=np.float)
        for t in range(1,self.t):
            for j in range(self.j):
                self.yub[t,j] = self.ub[t-1,j] + self.yub[t-1,j]
                self.ylb[t,j] = self.lb[t-1,j] + self.ylb[t-1,j]
        for j in range(self.j):
            self.yub[0,j] = 1
        self.produceMesh()
        tmp = {}
        cu = {}
        self.xi = np.zeros((self.t,self.j,self.a,self.b),dtype=np.float)
        for s in range(self.lenMonSec):
            #simulation
            for t in range(self.t):
                for j in range(self.j):
                    tmp[t,j]=np.sum(np.random.uniform(size=self.T)<self.prob[t*self.T:(t+1)*self.T,j])
            #calculate
            for j in range(self.j):
                cu[j] = 0
            for t in range(self.t):
                for j in range(self.j):
                    result = self.pieceWiseFunctionOnMesh(t,j,tmp[t,j],cu[j]+(t==0))
                    for (a,b,value) in result:
                        self.xi[t,j,a,b] += value
                    cu[j] += tmp[t,j]
        #mean
        for t in range(self.t):
            for j in range(self.j):
                for a in range(self.a):
                    for b in range(self.b):
                        self.xi[t,j,a,b] /= float(self.lenMonSec)        
        self.buildWh()
        #print self.xi
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
        return self.realDemand,tmpDemand
                
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
        