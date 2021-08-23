import numpy as np
from copy import copy, deepcopy
import networkx as nx
import random

from .solver import IK_Solver
from rcrilib.Helpers import IK_ParameterSet, IK_Parameter, ikdof, createLogger, IK_Math
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_FourFiveSolver")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class IK_FourFiveCycleSolver(IK_Solver):
    def __init__(self, G, ndof, consbonds, config, applycheck=True):
        self.G = G
        self.consbonds = consbonds
        self.config = config

        self.PS = IK_ParameterSet()
        for lb in consbonds:
            if len(lb.bond) == 2:
                self.PS += IK_Parameter(ikdof.CONTINUOUS, ikdof.DEPENDENT, copy(lb.side1), tc.SOLVER, atoms=copy(lb.bond))
        self.PS += IK_Parameter(ikdof.DISCRETE, ikdof.FREE, 0, tc.SOLVER)
        self.genconformers()

    def genconformers(self):
        N = len(self.G.nodes())
        self.nodes = list(self.G.nodes())
        self.conformers = [np.zeros(shape=(N,3)),np.zeros(shape=(N,3))]
        for i in range(N):
            self.conformers[0][i][:] = self.G.nodes[self.nodes[i]]['xyz'][:]
            self.conformers[1][i][:] = -self.G.nodes[self.nodes[i]]['xyz'][:]
        self.PS.getMyDiscreteParameter(tc.SOLVER).maxValue = 2

    @staticmethod
    def canapply(G, ndof, consbonds, config):
        if nx.number_of_nodes(G) < 4:
            return False

        if nx.number_of_nodes(G) != nx.number_of_edges(G): # NOT A SINGLE CYCLE
            return False

        for item in consbonds:
            if len(item.bond) == 2:
                idx = [item.side1[0], item.bond[0], item.bond[1], item.side1[1]]
                if abs(IK_Math.gettorsion([G.nodes[idx[0]]['xyz'],
                                    G.nodes[idx[1]]['xyz'],
                                    G.nodes[idx[2]]['xyz'],
                                    G.nodes[idx[3]]['xyz']])) > 44*deg2rad:
                    return False

        for edge in G.edges():
            if not G[edge[0]][edge[1]]['single']:
                leftlist = list(G.neighbors(edge[0]))
                leftlist.remove(edge[1])
                rightlist = list(G.neighbors(edge[1]))
                rightlist.remove(edge[0])
                vecs = [G.nodes[leftlist[0]]['xyz'], G.nodes[edge[0]]['xyz'],
                        G.nodes[edge[1]]['xyz'], G.nodes[rightlist[0]]['xyz']]
                if abs(IK_Math.gettorsion(vecs)) > 44*deg2rad:
                    return False
        return True

    def updatexyz(self, G):
        N = len(self.G.nodes())
        nodes = list(self.G.nodes())
        logger.info("Getting conformer #%d" % self.PS.getMyDiscreteParameter(tc.SOLVER).value)
        for i in range(N):
            G.nodes[nodes[i]]['xyz'][:] = self.conformers[self.PS.getMyDiscreteParameter(tc.SOLVER).value][i][:]

    def updateGraphxyz(self):
        N = len(self.G.nodes())
        nodes = list(self.G.nodes())
        logger.info("Getting conformer #%d" % self.PS.getMyDiscreteParameter(tc.SOLVER).value)
        for i in range(N):
            self.G.nodes[nodes[i]]['xyz'][:] = self.conformers[self.PS.getMyDiscreteParameter(tc.SOLVER).value][i][:]

    def getPS(self):
        return self.PS.getPS()

    def applyPS(self):
        # TODO Check if consbonds are complementary
        if self.PS.getMyDiscreteParameter(tc.SOLVER).value is None:
            self.PS.getMyDiscreteParameter(tc.SOLVER).value = random.randint(0, 1)
        confnum = self.PS.getMyDiscreteParameter(tc.SOLVER).value
        self.PS.getMyDiscreteParameter(tc.SOLVER).solutionCounter.recordState()

        for item in self.PS:
            if item.isDependent() and item.isContinuous():
                idx = [item.sides[0], item.atoms[0], item.atoms[1], item.sides[1]]
                for i in range(len(idx)):
                    idx[i] = self.nodes.index(idx[i])
                # TODO Check if torsion constraints are satisfied
        return IK_ParameterSet()
