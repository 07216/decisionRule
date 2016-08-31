# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 16:50:33 2016

@author: Zhan
"""


import Input
import sys
import decisionRule
import rALP
import time
import simulation
import CustomizeDemand
import gene
import numpy as np
#from guppy import hpy

firstCaseInReSolve=1
secondCaseInReSolve=2
reductionALP=0

demander = CustomizeDemand.CustomizeDemand(1,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))

first = 1

if first == 1:#Decision Rule Approch 
    
    decisionSolver = decisionRule.decisionRule()
    #   decisionSolver.echoInput()
    #decisionSolver.echoVal()
    decisionSolver.inputDemand(demander)
    start = time.clock()
    print "Input End",start
    decisionSolver.addVar()
    print "AddVar End",time.clock()
    decisionSolver.addOpt()
    print "AddOpt End",time.clock()
    decisionSolver.addConstr()
    end = time.clock()
    print "AddCons End",end
    #decisionSolver.writeMPS()
    decisionSolver.solve()
    #decisionSolver.echoOpt()    
    
    print end-start
    
    simulator = simulation.simulation(decisionSolver,demander)
    '''
    realDemand = {}
    for i in range(0,3):
        realDemand[i] = demander.sim()
    
    opt = gene.gene(decisionSolver,simulator)
    opt.initX()
    opt.setSTD(realDemand)
    opt.evolve()
    '''
    simulator.initX()
    simulator.initXX()
    #simulator.echoXX()
    print simulator.bookLimRun(np.round,1000)
    print simulator.bookLimRun(np.ceil,1000)
    print simulator.bookLimRun(np.floor,1000)
    print simulator.bookLimRun(simulator.atLeastOne,1000)
    print simulator.run(1000)

elif first == 0:#reduction of Approximate Linear Programming
    
    reduction = rALP.rALP(demander.reader)
    reduction.construct()
    reduction.addVar()
    reduction.addOpt()
    reduction.addConstr()
    reduction.solve()

elif first == 2:#Gene Algorithm
    
    decisionSolver = decisionRule.decisionRule()
    decisionSolver.inputDemand(demander)    
    decisionSolver.x = {}
    simulator = simulation.simulation(decisionSolver,demander)

    realDemand = {}
    for i in range(0,100):
        realDemand[i] = demander.sim()

    opt = gene.gene(decisionSolver,simulator)
    opt.zeroX()
    opt.setSTD(realDemand)
    opt.evolve()
    #simulator.run(1000)
    print simulator.bookLimRun(100,opt.X)
