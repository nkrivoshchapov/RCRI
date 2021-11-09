from copy import copy, deepcopy
import numpy as np
from numpy.linalg import inv
from numpy import cos, sin
import builtins

from rcrilib.Helpers import IK_Math, createLogger

logger = createLogger("IK_LinkingBond")

deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877
class IK_LinkingBond:
    def __init__(self, bond, side1, side2, ind_cycle_idx, dep_cycle_idx):
        # 1) TWO INTEGERS: ATOMS OF LINKING BOND, NUMERATION OF self.G
        # 2) TWO INTEGERS: ATOMS ON THE SIDE OF self.DG.nodes[node]
        # 3) TWO INTEGERS: ATOMS ON THE SIDE OF self.DG.nodes[nbnode]
        # node=side1 - dependent
        # nbnode=side2 - independent
        self.bond = copy(bond)
        self.side1 = copy(side1)
        self.side2 = copy(side2)
        self.orient = True  # TRUE - SAME ORIENTATION IN IK_LinkingBond and IK_TLCSolver, ELSE - FALSE
        self.ind_idx = ind_cycle_idx
        self.dep_idx = dep_cycle_idx
        self.done_init = False

    def __eq__(self, other):
        if isinstance(other, IK_LinkingBond):
            return self.bond == other.bond and self.side1 == other.side1 and self.side2 == other.side2
        elif isinstance(other, list): # Duplicate of hasBond
            return other[0] in self.bond and other[1] in self.bond
        else:
            raise Exception("Invalid comparison")

    def hasBond(self, atompair):
        if isinstance(atompair, list):
            return atompair[0] in self.bond and atompair[1] in self.bond
        else:
            raise Exception("Invalid comparison")

    def __len__(self):
        return len(self.bond)

    def setParam(self, param):
        self.param = param
        if not param.isDependent() and not param.isFixed():
            raise Exception('Trying to establish dependence of independent parameter!')

    def readDependence(self, G):
        self.done_init = True
        if len(self.bond) == 1:
            #SIDE 1
            at0 = G.nodes[self.side1[0]]['xyz']
            at1 = G.nodes[self.bond[0]]['xyz']
            at2 = G.nodes[self.side1[1]]['xyz']
            xv1 = at2 - at1
            yv1 = at0 - at1
            zv1 = IK_Math.gs_rand(xv1, yv1)
            frame1 = np.array([xv1, yv1, zv1]).transpose()

            # SIDE 2
            at0_2 = G.nodes[self.side2[0]]['xyz']
            at2_2 = G.nodes[self.side2[1]]['xyz']
            xv2 = at2_2 - at1
            yv2 = at0_2 - at1
            zv2 = IK_Math.gs_rand(xv2, yv2)
            frame2 = np.array([xv2, yv2, zv2]).transpose()
            self.transit_mat = inv(frame1) @ frame2
        elif len(self.bond) == 2:
            # SIDE 1
            at0 = G.nodes[self.side1[0]]['xyz']
            at1 = G.nodes[self.bond[0]]['xyz']
            at2 = G.nodes[self.bond[1]]['xyz']
            at3 = G.nodes[self.side1[1]]['xyz']
            tor1 = IK_Math.gettorsion([at0, at1, at2, at3])

            # SIDE 2
            at0_2 = G.nodes[self.side2[0]]['xyz']
            at3_2 = G.nodes[self.side2[1]]['xyz']
            tor2 = IK_Math.gettorsion([at0_2, at1, at2, at3_2])

            self.tordiff = tor1 - tor2  # self.tordiff MUST BE ADDED TO INDEPENDENT TORSION!
            self.twosidediff = -IK_Math.gettorsion([at0, at1, at2, at0_2])
            self.transit_mat = np.array([[1, 0, 0],
                                         [0, cos(self.twosidediff), -sin(self.twosidediff)],
                                         [0, sin(self.twosidediff), cos(self.twosidediff)]])
            """
            EQUIVALENT WAY TO FIND self.transit_mat 
            xv1 = at2 - at1
            xv1 /= norm(xv1)
            zv1 = np.cross(xv1, at0 - at1)
            zv1 /= norm(zv1)
            yv1 = np.cross(zv1, xv1)
            frame1 = np.array([xv1, yv1, zv1]).transpose()
            xv2 = at2 - at1
            xv2 /= norm(xv2)
            zv2 = np.cross(xv2, at0_2 - at1)
            zv2 /= norm(zv2)
            yv2 = np.cross(zv2, xv2)
            frame2 = np.array([xv2, yv2, zv2]).transpose()
            self.transit_mat = inv(frame1) @ frame2
            """

    def resolveDependence(self, ind_G):
        if len(self.bond) == 2:
            at0 = ind_G.nodes[self.side2[0]]['xyz']
            at1 = ind_G.nodes[self.bond[0]]['xyz']
            at2 = ind_G.nodes[self.bond[1]]['xyz']
            at3 = ind_G.nodes[self.side2[1]]['xyz']
            logger.debug("My atom numbers side1(dependent) = " +
                  repr([self.side1[0], self.bond[0], self.bond[1], self.side1[1]]))
            logger.debug("My atom numbers side2(independent) = " +
                  repr([self.side2[0], self.bond[0], self.bond[1], self.side2[1]]))
            if self.orient:
                logger.debug("Positive sign:" + repr([self.side2[0], self.bond[0],
                                                             self.bond[1], self.side2[1]]))
                tor1 = (self.tordiff + IK_Math.gettorsion([at0, at1, at2, at3]))  # QUESTIONABLE SIGN
            else:
                tor1 = (IK_Math.gettorsion([at0, at1, at2, at3]) + self.tordiff)  # QUESTIONABLE SIGN
                logger.debug("Negative sign:" + repr([self.side2[0], self.bond[0],
                                                             self.bond[1], self.side2[1]]))
            if tor1 > np.pi:
                tor1 -= 2 * np.pi
            if tor1 < -np.pi:
                tor1 += 2 * np.pi
            logger.debug("My dep " + repr([self.side2[0], self.bond[0], self.bond[1], self.side2[1]]) +
                  " torsion = " + repr(tor1 * rad2deg))
            self.param.setValue(tor1)

    def checkOrientation(self):
        bond = [self.param.atoms[0], self.param.atoms[1]]
        logger.debug("The parameter bond is" + repr(bond))
        logger.debug("My bond is" + repr(self.bond))
        if bond[0] == self.bond[1] and bond[1] == self.bond[0]:
            logger.debug("Negative orientation")
            self.orient = False
        elif not (bond[1] == self.bond[1] and bond[0] == self.bond[0]):
            raise Exception("Parameter and Linking bond don't match")

    def connect_cycles(self, G_dep, G_ind, G_glob):
        if len(self.bond) == 1:
            at0 = G_dep.nodes[self.side1[0]]['xyz']
            at1 = G_dep.nodes[self.bond[0]]['xyz']
            at2 = G_dep.nodes[self.side1[1]]['xyz']
            av = at2 - at1
            bv = at0 - at1
            cv = IK_Math.gs_rand(av, bv)

            framemat = np.zeros((4, 4))
            framemat[:3, :3] = np.array([av, bv, cv]).transpose() @ self.transit_mat
            framemat[3][3] = 1
            framemat[:3, 3] = at1
            baseat = [self.bond[0], self.side2[1], self.side2[0]]
        elif len(self.bond) == 2:
            at0 = G_dep.nodes[self.side1[0]]['xyz']
            at1 = G_dep.nodes[self.bond[0]]['xyz']
            at2 = G_dep.nodes[self.bond[1]]['xyz']

            av = at2 - at1
            bv = at0 - at1
            cv = IK_Math.gs_rand(av, bv)
            
            framemat = np.zeros((4, 4))
            framemat[:3, :3] = np.array([av, bv, cv]).transpose() @ self.transit_mat
            framemat[3][3] = 1
            framemat[:3, 3] = at1
            baseat = [self.bond[0], self.bond[1], self.side2[0]]
        else:
            raise Exception("How many linking atoms?")

        center = deepcopy(G_ind.nodes[baseat[0]]['xyz'])
        av = G_ind.nodes[baseat[1]]['xyz'] - center
        bv = G_ind.nodes[baseat[2]]['xyz'] - center
        cv = IK_Math.gs_rand(av, bv)

        for node in list(G_ind.nodes()):
            coord = framemat @ np.array([(G_ind.nodes[node]['xyz'] - center) @ av,
                                         (G_ind.nodes[node]['xyz'] - center) @ bv,
                                         (G_ind.nodes[node]['xyz'] - center) @ cv,
                                         1])
            G_glob.nodes[node]['xyz'] = np.array([coord[0], coord[1], coord[2]])

        # xyzlines = [str(G_ind.number_of_nodes()), ""]
        # for node in list(G_ind.nodes()):
        #     xyzlines.append("%3s%10.4f%10.4f%10.4f%3s" % ("C", G_ind.nodes[node]['xyz'][0],
        #                                                   G_ind.nodes[node]['xyz'][1],
        #                                                   G_ind.nodes[node]['xyz'][2], ""))
        # wfile = open("indepframe2.xyz", "w")
        # wfile.write("\n".join(xyzlines))
        # wfile.close()

        # if builtins.loglinecount == 66:
        #     print("HERE")

    def isFake(self):
        return False
