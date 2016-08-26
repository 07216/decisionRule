# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 16:50:33 2016

@author: Zhan
"""


import Input
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

demander = CustomizeDemand.CustomizeDemand(2)

first = 1

if first == 1:#Decision Rule Approch 
    
    decisionSolver = decisionRule.decisionRule()
    #   decisionSolver.echoInput()
    #decisionSolver.echoVal()
    decisionSolver.inputDemand(demander)
    decisionSolver.addVar()
    decisionSolver.addOpt()
    decisionSolver.addConstr()
    #decisionSolver.writeMPS()
    decisionSolver.solve()
    #decisionSolver.echoOpt()    
    
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
    simulator.initXX(np.floor,0)
    simulator.echoXX()
    print simulator.run(100)
    simulator.initXX(np.floor,0)
    print simulator.bookLimRun(1000)
    simulator.initXX(np.ceil,0)
    print simulator.bookLimRun(1000)
    simulator.initXX(simulator.atLeastOne,0)
    print simulator.bookLimRun(1000)
    simulator.initXX(np.round,0)
    print simulator.bookLimRun(1000)
    simulator.initXX(np.floor,1)
    print simulator.bookLimRun(1000)
    simulator.initXX(np.ceil,1)
    print simulator.bookLimRun(1000)
    simulator.initXX(simulator.atLeastOne,1)
    print simulator.bookLimRun(1000)
    simulator.initXX(np.round,1)
    print simulator.bookLimRun(1000)

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
