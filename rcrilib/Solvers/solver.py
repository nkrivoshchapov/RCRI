import numpy as np
from copy import copy

from rcrilib.Helpers import IK_ParameterSet, IK_Parameter, ikdof, cause, createLogger, IK_Math
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_IdentitySolver")

class IK_Solver:
    @staticmethod
    def add_constraints(G):
        for edge in list(G.edges()):
            G[edge[0]][edge[1]]['type'] = False

class IK_IdentitySolver(IK_Solver):
    def __init__(self, G, ndof, consbonds, config, applycheck=True):
        self.G = G
        self.consbonds = consbonds
        self.config = config

        self.PS = IK_ParameterSet()
        for lb in consbonds:
            if len(lb.bond) == 2:
                self.PS += IK_Parameter(ikdof.CONTINUOUS, ikdof.DEPENDENT, copy(lb.side1),
                                        tc.SOLVER, atoms=copy(lb.bond))

        N = len(self.G.nodes())
        self.nodes = list(self.G.nodes())
        self.conformer = np.zeros(shape=(N, 3))
        for i in range(N):
            self.conformer[i][:] = self.G.nodes[self.nodes[i]]['xyz'][:]

    def updatexyz(self, G):
        for at in list(self.G.nodes):
            G.nodes[at]['xyz'] = np.array([self.G.nodes[at]['xyz'][0],
                                           self.G.nodes[at]['xyz'][1],
                                           self.G.nodes[at]['xyz'][2]])

    def updateGraphxyz(self):
        pass # BECAUSE ONLY THIS SOLVER PRESERVES CONFORMATION. updateGraphxyz MUST BE IMPLEMENTED FOR ALL OTHER SOLVERS

    def getPS(self):
        return self.PS.getPS()

    def applyPS(self):
        self.PS.freeze()
        for item in self.PS:
            if item.isDependent() and item.isContinuous():
                idx = [item.sides[0], item.atoms[0], item.atoms[1], item.sides[1]]
                for i in range(len(idx)):
                    idx[i] = self.nodes.index(idx[i])
                logger.debug("Required = %f; Actual = %f" % (item.value, IK_Math.gettorsion([self.conformer[idx[0]],
                                                                                     self.conformer[idx[1]],
                                                                                     self.conformer[idx[2]],
                                                                                     self.conformer[idx[3]]])))
                if abs(item.value - IK_Math.gettorsion([self.conformer[idx[0]],
                                                     self.conformer[idx[1]],
                                                     self.conformer[idx[2]],
                                                     self.conformer[idx[3]]]) > 0.001):
                    self.PS.cause = cause.zerosolutions
                    return self.PS.getPS(excludeFixed=True, errorcause=cause.zerosolutions_ddof)

        return IK_ParameterSet()
