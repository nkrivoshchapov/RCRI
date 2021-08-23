import numpy as np
from copy import deepcopy
from numpy.linalg import norm
from enum import Enum
from scipy import spatial

from .math import IK_Math
from .logging import createLogger

logger = createLogger("Validation")

class constype(Enum):
    LENGTH = 1
    VANGLE = 2
    TANGLE = 3
    POLYHEDRA = 4

deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class IK_GC_BondLength:
    def __init__(self, G, idx):
        self.type = constype.LENGTH
        self.idx = deepcopy(idx)
        # idx = [Gidx of first atom,Gidx of second atom]
        self.value = norm(G.nodes()[idx[0]]['xyz']-G.nodes()[idx[1]]['xyz'])

    def validate(self, G, attr="xyz"):
        if not abs(norm(G.nodes()[self.idx[0]][attr]-G.nodes()[self.idx[1]][attr])-self.value) < 0.001:
            logger.error("Constraint on bond %s is NOT satisfied: %f instead of %f" % (repr(self.idx),
                                        norm(G.nodes()[self.idx[0]][attr]-G.nodes()[self.idx[1]][attr]),self.value))
            return False
        else:
            logger.debug("Constraint on bond %s is satisfied: %f instead of %f" % (repr(self.idx),
                                                                                norm(G.nodes()[self.idx[0]][attr] -
                                                                                     G.nodes()[self.idx[1]][attr]),
                                                                                self.value))
            return True

class IK_GC_ValAngle:
    def __init__(self, G, idx):
        self.type = constype.VANGLE
        self.idx = deepcopy(idx)
        # ind = [Gidx of side atom 1, Gidx of middle atom, Gidx of side atom 2]
        self.value = IK_Math.getvalangle([G.nodes()[idx[0]]['xyz'],G.nodes()[idx[1]]['xyz'],G.nodes()[idx[2]]['xyz']])

    def validate(self, G, attr="xyz"):
        if not abs(IK_Math.getvalangle([G.nodes()[self.idx[0]][attr],G.nodes()[self.idx[1]][attr],G.nodes()[self.idx[2]][attr]])-self.value) < 0.001:
            logger.error("Constraint on vangle %s is NOT satisfied: %f instead of %f" % (repr(self.idx),
                IK_Math.getvalangle([G.nodes()[self.idx[0]][attr],G.nodes()[self.idx[1]][attr],G.nodes()[self.idx[2]][attr]]),
                self.value))
            return False
        else:
            logger.debug("Constraint on vangle %s is satisfied: %f instead of %f" % (repr(self.idx),
              IK_Math.getvalangle([G.nodes()[self.idx[0]][attr], G.nodes()[self.idx[1]][attr], G.nodes()[self.idx[2]][attr]]),
              self.value))
            return True

class IK_GC_Polyhedra:
    def __init__(self, G, idx):
        self.type = constype.POLYHEDRA
        self.idx = deepcopy(idx)
        # idx = [[Gidx of side atom 1, Gidx of middle atom, Gidx of side atom 2],[other atoms]]
        prevvec = G.nodes()[idx[0][0]]['xyz']
        curvec = G.nodes()[idx[0][1]]['xyz']
        nextvec = G.nodes()[idx[0][2]]['xyz']

        av = nextvec - curvec
        bv = prevvec - curvec
        cv = IK_Math.gs_rand(av, bv)

        self.value = []
        for item in idx[1]:
            temp = G.nodes()[item]['xyz'] - curvec
            self.value.append(np.array([np.dot(temp, av), np.dot(temp, bv), np.dot(temp, cv)]))

    def validate(self, G, attr="xyz"):
        prevvec = G.nodes()[self.idx[0][0]][attr]
        curvec = G.nodes()[self.idx[0][1]][attr]
        nextvec = G.nodes()[self.idx[0][2]][attr]

        av = nextvec - curvec
        bv = prevvec - curvec
        cv = IK_Math.gs_rand(av, bv)

        res = True
        for i in range(len(self.value)):

            temp = G.nodes()[self.idx[1][i]][attr] - curvec
            temp = np.array([np.dot(temp,av),np.dot(temp,bv),np.dot(temp,cv)])
            if not norm(self.value[i] - temp)<0.005:
                logger.error("Constraint on position of %d in polyhedra %d is NOT satisfied: %s instead of %s "
                      "(misplacement = %f)" % (self.idx[1][i], self.idx[0][1],
                                               repr(self.value[i]), repr(temp), norm(self.value[i] - temp)))
                res = False
            else:
                logger.debug("Constraint on position of %d in polyhedra %d is satisfied: %s instead of %s "
                      "(misplacement = %f)" % (self.idx[1][i], self.idx[0][1],
                                               repr(self.value[i]), repr(temp), norm(self.value[i] - temp)))
        return res

class IK_GC_Torsion:
    def __init__ (self, G, idx, param):
        self.type = constype.TANGLE
        self.idx = deepcopy(idx)
        self.value = 0
        self.param = param
        # ind = [Gidx of side atom 1, Gidx of middle atom 1, Gidx of middle atom 2,  Gidx of side atom 2]

    def validate(self, G, attr="xyz"):
        curvalue = IK_Math.gettorsion([G.nodes()[self.idx[0]][attr],G.nodes()[self.idx[1]][attr],
                                       G.nodes()[self.idx[2]][attr],G.nodes()[self.idx[3]][attr]])
        if not abs(curvalue-self.param.value) < 0.001:
            logger.error("Constraint on tangle %s is NOT satisfied: %f instead of %f" % (repr(self.idx),
                curvalue * rad2deg, self.param.value * rad2deg))
            return False
        else:
            logger.debug("Constraint on tangle %s is satisfied: %f instead of %f" % (repr(self.idx),
                curvalue * rad2deg, self.param.value * rad2deg))
            return True

class IK_GeomValidator:
    def __init__(self, outerG, PS):
        self.G = outerG
        self.constraints = []
        self.initialize(PS)

    def initialize(self, PS):
        #BONDS
        for bond in list(self.G.edges):
            self.constraints.append(IK_GC_BondLength(self.G,bond))

        # #VAL ANGLES AND POLYHEDRA
        # #TODO get rid of linearity problem
        # for node in list(self.G.nodes()):
        #     nb = list(self.G.neighbors(node))
        #     if len(nb) >= 2:
        #         self.constraints.append(IK_GC_ValAngle(self.G,[nb[0],node,nb[1]]))
        #         if len(nb) > 2:
        #             othernodes = []
        #             for i in range(2,len(nb)):
        #                 othernodes.append(nb[i])
        #             self.constraints.append(IK_GC_Polyhedra(self.G, [[nb[0], node, nb[1]],othernodes]))
        # TORSIONS
        # for item in PS:
        #     if item.isContinuous() and len(item.atoms) == 2:
        #         self.constraints.append(IK_GC_Torsion(self.G, [item.sides[0], item.atoms[0], item.atoms[1], item.sides[1]], item))

    def validate(self, G, attr="xyz"):
        passed = True
        for cons in self.constraints:
            ok = cons.validate(G, attr=attr)
            if not ok:
                passed = False
                raise Exception("Bond length check not passed")
        return passed

def geomcheckTree_array(atoms_xyz, bonds, mindist):
    tree = spatial.KDTree(atoms_xyz)
    for i, nbs in enumerate(tree.query_ball_point([atoms_xyz], mindist)[0]):
        for nb in nbs:
            if nb != i and [i, nb] not in bonds and [nb, i] not in bonds:
                logger.info("Overlap detected")
                return False
    return True

def geomcheckTree_graph(graph, mindist, attr='xyz'):
    atoms_xyz = np.zeros((graph.number_of_nodes(), 3))
    names = []
    for i, node in enumerate(list(graph.nodes())):
        names.append(node)
        atoms_xyz[i] = deepcopy(graph.nodes[node][attr])

    tree = spatial.KDTree(atoms_xyz)
    for i, nbs in enumerate(tree.query_ball_point([atoms_xyz], mindist)[0]):
        for nb in nbs:
            if nb != i and names[nb] not in graph.neighbors(names[i]):
                logger.info("Overlap detected")
                return False
    return True

def geomcheckTree_cycle(cycle, mindist, attr='xyz'):
    atoms_xyz = np.zeros((cycle.ord, 3))
    names = []
    for i, node in enumerate(cycle.items):
        names.append(node)
        atoms_xyz[i] = deepcopy(cycle.items[i][attr])
    tree = spatial.KDTree(atoms_xyz)
    for i, nbs in enumerate(tree.query_ball_point([atoms_xyz], mindist)[0]):
        for nb in nbs:
            if nb != i and nb != (cycle.items[i]-1)['num'] and nb != (cycle.items[i]+1)['num']:
                logger.info("Overlap detected")
                return False
    return True

def geomcheckTree_cycpart(graph, extgraph, mindist, attr='xyz'):
    nnodes = graph.number_of_nodes()
    for cnode in list(graph.nodes()):
        nnodes += len(graph.nodes[cnode]['outerNB'])
    atoms_xyz = np.zeros((nnodes, 3))
    names = []
    j = 0
    for cnode in list(graph.nodes()):
        names.append(cnode)
        atoms_xyz[j] = deepcopy(graph.nodes[cnode][attr])
        j += 1
    for cnode_idx, cnode in enumerate(list(graph.nodes())):
        for node_idx, node in enumerate(graph.nodes[cnode]['outerNB']):
            names.append(node)
            atoms_xyz[j] = deepcopy(graph.nodes[cnode]['outerXYZ'][node_idx])
            j += 1

    tree = spatial.KDTree(atoms_xyz)
    for i, nbs in enumerate(tree.query_ball_point([deepcopy(atoms_xyz)], mindist)[0]):
        for nb in nbs:
            if nb != i and names[nb] not in extgraph.neighbors(names[i]):
                logger.info("Overlap detected")
                return False
    return True
