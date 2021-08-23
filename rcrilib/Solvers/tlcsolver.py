from ctypes import *
from copy import copy,deepcopy
import networkx as nx
import numpy as np
import os, random

from .solver import IK_Solver
from rcrilib.Helpers import Cngroup, GroupItem, itemtype, IK_ParameterSet, IK_Parameter, ikdof, getfloat, getbool, cause,\
    geomcheckTree_cycle, createLogger, IK_Math
from .rigidfrag import IK_RigidFrag
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_TLCSolver")

deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class TLC_xyz_arr(Structure):
    _fields_ = [('nconf', c_int), ('arr', ((c_double*3)*9)*16)]

class TLC_str_data(Structure):
    _fields_ = [('leng', (c_double*9)), ('vang', (c_double*9)), ('tang', (c_double*3))]

class IK_TLCsolver(IK_Solver):
    def __init__(self, G, ndof, consbonds, config):
        self.G = G
        self.config = config
        self.gengroups(consbonds)
        self.initbonds()
        self.initrigidfrags()

        basedir = os.path.abspath(os.path.dirname(__file__))
        TLC_libfile = os.path.join(basedir, 'libtlc.so')
        TLC_dll = CDLL(TLC_libfile)
        self.SolveTLC = TLC_dll.genconf
        self.SolveTLC.argtypes = [POINTER(TLC_str_data), POINTER(TLC_xyz_arr)]

    @staticmethod
    def walk_cycle(G):
        walk = [list(G.nodes())[0]]
        curnode = list(G.neighbors(walk[0]))[0]
        i = 1
        while curnode != walk[0]:
            walk.append(curnode)
            if list(G.neighbors(curnode))[0] == walk[i-1]:
                curnode = list(G.neighbors(curnode))[1]
            else:
                curnode = list(G.neighbors(curnode))[0]
            i += 1
        return walk

    def gengroups(self, consbonds):
        seq_atoms = IK_TLCsolver.walk_cycle(self.G)
        logger.info("TLC applied for cycle of %d atoms" % len(seq_atoms))
        atoms = []
        prevatom = -1
        curatom = seq_atoms[0]
        self.Catoms = Cngroup(len(seq_atoms), itype=itemtype.ATOM)
        self.Cbonds = Cngroup(len(seq_atoms), itype=itemtype.BOND)  # i-th bond is between i and i+1 atoms
        i = 0
        while (not curatom == seq_atoms[0]) or (prevatom == -1):
            atoms.append(curatom)
            if prevatom == -1:
                temp = list(self.G.neighbors(curatom))[0]
            else:
                singleatom = list(self.G.neighbors(curatom))
                singleatom.remove(prevatom)
                temp = singleatom[0]
            prevatom = curatom
            curatom = temp

            mynext = i + 1
            if mynext == len(seq_atoms):
                mynext = 0
            # current bond is between i(=prevatom) and mynext(=curatom)
            self.Catoms[i]['G_idx'] = prevatom
            self.Catoms[i]['bonds'] = [self.Cbonds[i - 1], self.Cbonds[i]]
            self.Cbonds[i]['atoms'] = [self.Catoms[i], self.Catoms[i + 1]]

            # BONDS ARE EITHER ROTATABLE OR NONROTATABLE. NONROTATABLE MAY HAVE CONSTANT TORSION (double bonds, etc.)
            # OR TORSION WILL BE DETERMINED BY CONFORMATION OF ANOTHER CYCLE

            if self.G[curatom][prevatom]['type']:
                if [curatom, prevatom] in consbonds or [prevatom, curatom] in consbonds:
                    self.Cbonds[i]['free'] = False
                    self.Cbonds[i]['dep'] = True
                else:
                    self.Cbonds[i]['dep'] = False
                    self.Cbonds[i]['free'] = True
            else:
                self.Cbonds[i]['free'] = False
                self.Cbonds[i]['dep'] = False
            i = mynext

    def initbonds(self):
        startbond = self.Cbonds.bond_getstart()
        mybond = startbond
        atompos = []
        while True:
            if mybond['free'] and (mybond + 1)['free']:
                myatom = mybond.bond_getcommonatom(mybond + 1)
                atompos.append(myatom)
                mybond += 2
            else:
                mybond += 1

            if len(atompos) == 3 or mybond == startbond:
                break
        logger.debug("STARTING CHOICE OF THE FRAME:")
        for i in range(len(atompos)):
            logger.debug("%s : %s" % (repr(atompos[i]), repr(atompos[i]['bonds'])))
        logger.debug("With metric = " + repr(self.getmetric(atompos[0], atompos[1], atompos[2])))

        curmetric = -1
        pa = self.Catoms.atoms_find_proper_atoms()
        for i in range(len(pa)):
            for j in range(i + 1, len(pa)):
                for k in range(j + 1, len(pa)):
                    if self.check_if_sensible(pa[i], pa[j], pa[k]):
                        metric = self.getmetric(pa[i], pa[j], pa[k])
                        if metric > curmetric:
                            atompos[0] = pa[i]
                            atompos[1] = pa[j]
                            atompos[2] = pa[k]
                            curmetric = metric
        self.axatoms = copy(atompos)
        logger.debug("BEST CHOICE OF THE FRAME:")
        for i in range(len(atompos)):
            logger.debug("%s : %s" % (repr(atompos[i]), repr(atompos[i]['bonds'])))
        logger.debug("With metric = " + repr(self.getmetric(atompos[0], atompos[1], atompos[2])))

    def getmetric(self, a, b, c):
        return min(self.Catoms.atoms_dist(a, b),
                   self.Catoms.atoms_dist(b, c),
                   self.Catoms.atoms_dist(a, c))

    def check_if_sensible(self, a, b, c):
        return self.getmetric(a, b, c) > 0

    @staticmethod
    def canapply(G, ndof, consbonds):
        if ndof < 0:
            return False

        if nx.number_of_nodes(G) != nx.number_of_edges(G): # NOT A SINGLE CYCLE
            return False
        atoms = []
        bondstring = ""
        prevatom = -1
        mcb = IK_TLCsolver.walk_cycle(G)
        curatom = mcb[0]
        while (not curatom == mcb[0]) or (prevatom == -1):
            atoms.append(curatom)
            if prevatom == -1:
                temp = list(G.neighbors(curatom))[0]
            else:
                singleatom = list(G.neighbors(curatom))
                singleatom.remove(prevatom)
                temp = singleatom[0]
            prevatom = curatom
            curatom = temp
        i = 0
        while True:
            ed1 = atoms[i]
            if i == len(atoms) - 1:
                ed2 = atoms[0]
            else:
                ed2 = atoms[i + 1]

            if G[ed1][ed2]['type'] and [ed1, ed2] not in consbonds and [ed2, ed1] not in consbonds:
                bondstring += "T"
            else:
                bondstring += "F"
            i += 1
            if i == len(atoms):
                break

        if "F" in bondstring:
            prev_letter = "F"
            count = 0
            for i in list(range(bondstring.index("F") + 1, len(bondstring))) + list(range(0, bondstring.index("F"))):
                if prev_letter == "T" and bondstring[i] == "T":
                    prev_letter = "F"
                    count += 1
                    if count == 3:
                        break
                else:
                    prev_letter = bondstring[i]
            if count >= 3:
                return True
            else:
                return False
        else:
            if len(bondstring) >= 6:
                return True
            else:
                logger.error("DOF Error")
                return False

    def initrigidfrags(self):
        for atom in self.Catoms.items:
            atom['xyz'] = deepcopy(self.G.nodes[atom['G_idx']]['xyz'])
        self.RF = []
        for pair in [[self.axatoms[0], self.axatoms[1]], [self.axatoms[1], self.axatoms[2]],
                     [self.axatoms[2], self.axatoms[0]]]:
            step = self.Catoms.atoms_getstep(pair, self.axatoms)
            if step == -1:
                temp = pair[0]
                pair[0] = pair[1]
                pair[0] = temp
            newfrag = IK_RigidFrag(pair, self.config)
            self.RF.append(newfrag)
            logger.debug("New fragment: %s" % (repr(newfrag)))

        self.PS = IK_ParameterSet()
        for item in self.RF:
            self.PS += item.getPS()
        self.PS += IK_Parameter(ikdof.DISCRETE, ikdof.FREE, 0, tc.SOLVER)
        logger.debug("PS representation in solver: " + repr(self.PS))

        # Precalculate angles
        self.RF[0].eatoms[0].atom_getangle()
        self.RF[1].eatoms[0].atom_getangle()
        self.RF[2].eatoms[0].atom_getangle()

    def getPS(self):
        return self.PS.getPS()

    def writecyclegeom(self):
        sat = self.axatoms[0]
        at = self.axatoms[0]
        while True:
            logger.debug("%d %14.6f %14.6f %14.6f" % (at['G_idx'], at['xyz'][0], at['xyz'][1], at['xyz'][2]))
            at += 1
            if at == sat:
                break

    def checkgeom(self):
        logger.info("Running geometry check...")
        passed = True
        for atom in self.Catoms.items:
            logger.debug("checking atom "+repr(atom))
            if not abs(atom.atom_getangle(forcecalc=True)-atom.atom_getangle(forcecalc=False)) < 0.001:
                passed = False
                logger.error("VAngle %s constraint is not satisfied - %f instead of %f" % (repr(atom),
                                            atom.atom_getangle(forcecalc=True) * rad2deg, atom.atom_getangle(forcecalc=False) * rad2deg))

        for bond in self.Cbonds.items:
            if not abs(bond.bond_getlength(forcecalc=True) - bond.bond_getlength(forcecalc=False)) < 0.001:
                passed = False
                logger.error("Bondlength %s constraint is not satisfied - %f instead of %f" % (repr(bond),
                                            bond.bond_getlength(forcecalc=True),bond.bond_getlength(forcecalc=False)))
        if not passed:
            raise Exception("[TLCSolver] Not a solution")
        else:
            logger.info("Test passed")

    def applyPS(self):
        if getbool("SmartAPS", "IK_TLCSolver", self.config) and not self.PS.isChanged():
            logger.info("Skipping TLC")
            if self.PS.cause == cause.success:
                return IK_ParameterSet()
            else:
                return self.PS.getPS(excludeFixed=True, errorcause=self.PS.cause)

        bonds = [self.RF[0].ebonds[0].bond_getlength(), None, self.RF[0].ebonds[1].bond_getlength(),
                 self.RF[1].ebonds[0].bond_getlength(), None, self.RF[1].ebonds[1].bond_getlength(),
                 self.RF[2].ebonds[0].bond_getlength(), None, self.RF[2].ebonds[1].bond_getlength()]
        vangles = [self.RF[0].eatoms[0].atom_getangle() * rad2deg, None, None,
                   self.RF[1].eatoms[0].atom_getangle() * rad2deg, None, None,
                   self.RF[2].eatoms[0].atom_getangle() * rad2deg, None, None]

        tangles = []
        for frag in self.RF:
            frag.applyPS()
            tangles.append(-frag.get_tangle())

        bonds[1] = self.RF[0].get_length()
        bonds[4] = self.RF[1].get_length()
        bonds[7] = self.RF[2].get_length()
        vangles[1], vangles[2] = self.RF[0].get_vangles()
        vangles[4], vangles[5] = self.RF[1].get_vangles()
        vangles[7], vangles[8] = self.RF[2].get_vangles()
        logger.debug("My bonds = %s\nMy vangles = %s\nMy tangles = %s" % (repr(bonds), repr(vangles), repr(tangles)))

        self.PS.freeze()

        xyzout = TLC_xyz_arr()
        inpdata = TLC_str_data()
        for i in range(9):
            inpdata.leng[i] = bonds[i]
            inpdata.vang[i] = vangles[i]
        for i in range(3):
            inpdata.tang[i] = tangles[i]
        self.SolveTLC(byref(inpdata), byref(xyzout))
        logger.info("We have %d solutions!" % (xyzout.nconf))
        confs = []
        for c in range(xyzout.nconf):
            conf = []
            for i in range(9):
                subj = []
                for j in range(3):
                    subj.append(xyzout.arr[c][i][j])
                conf.append(copy(subj))
            confs.append(deepcopy(conf))
        if len(confs) == 0:
            self.PS.cause = cause.zerosolutions
            if self.PS.getMyDiscreteParameter(tc.SOLVER).value is not None:
                raise Exception("Nonempty set of TLC solutions was expected")
            return self.PS.getPS(excludeFixed=True, errorcause=cause.zerosolutions)

        self.PS.getMyDiscreteParameter(tc.SOLVER).maxValue = len(confs)
        if self.PS.getMyDiscreteParameter(tc.SOLVER).value is None:
            self.PS.getMyDiscreteParameter(tc.SOLVER).value = random.randint(0, len(confs) - 1)
            logger.info("Randomized solution " + repr(self.PS.getMyDiscreteParameter(tc.SOLVER).value))
        self.PS.getMyDiscreteParameter(tc.SOLVER).solutionCounter.recordState()
        if not self.quickcheck(confs[self.PS.getMyDiscreteParameter(tc.SOLVER).value]):
            self.PS.cause = cause.zerosolutions
            return self.PS.getPS(excludeFixed=True, errorcause=cause.zerosolutions_ddof)
        self.setgeom(confs[self.PS.getMyDiscreteParameter(tc.SOLVER).value])

        logger.info("Getting solution #" + repr(self.PS.getMyDiscreteParameter(tc.SOLVER).value))
        for frag in self.RF:
            frag.setgeom()

        if getbool("DoValidation", "IK_TLCSolver", self.config):
            self.checkgeom()
        self.PS.cause = cause.success
        return IK_ParameterSet()

    def quickcheck(self, conf):
        ang1 = [np.array(conf[8]), np.array(conf[0]), np.array(conf[1])]
        ang2 = [np.array(conf[2]), np.array(conf[3]), np.array(conf[4])]
        ang3 = [np.array(conf[5]), np.array(conf[6]), np.array(conf[7])]

        if abs(IK_Math.getvalangle(ang1) - self.RF[0].eatoms[0].atom_getangle()) < 0.001 and \
               abs(IK_Math.getvalangle(ang2) - self.RF[1].eatoms[0].atom_getangle()) < 0.001 and \
               abs(IK_Math.getvalangle(ang3) - self.RF[2].eatoms[0].atom_getangle()) < 0.001:
            return True
        else:
            logger.warning("TLC failure detected")
            return False

    def setgeom(self, conf):
        self.axatoms[0]['xyz'] = np.array(conf[0])
        self.axatoms[1]['xyz'] = np.array(conf[3])
        self.axatoms[2]['xyz'] = np.array(conf[6])

        (self.axatoms[0] + 1)['xyz'] = np.array(conf[1])
        (self.axatoms[1] + 1)['xyz'] = np.array(conf[4])
        (self.axatoms[2] + 1)['xyz'] = np.array(conf[7])

        (self.axatoms[1] - 1)['xyz'] = np.array(conf[2])
        (self.axatoms[2] - 1)['xyz'] = np.array(conf[5])
        (self.axatoms[0] - 1)['xyz'] = np.array(conf[8])

    def updatexyz(self, G):
        for atom in self.Catoms.items:
            G.nodes[atom['G_idx']]['xyz'][:] = atom['xyz'][:]

    def updateGraphxyz(self):
        for atom in self.Catoms.items:
            self.G.nodes[atom['G_idx']]['xyz'][:] = atom['xyz'][:]

    def checkoverlap(self):
        if getbool("RunGeomCheck", "IK_TLCSolver", self.config) and \
                not geomcheckTree_cycle(self.Catoms, getfloat("MinDistance", "IK_TLCSolver", self.config)):
            return self.PS.getPS(excludeFixed=True, errorcause=cause.geomoverlap)

    def perturb_geometry(self):
        if getbool("AllowAnglePerturbation", "IK_TLCSolver", self.config):
            self.Catoms.perturb_angles(getfloat("AngleSDeviation", "IK_TLCSolver", self.config))

        if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config):
            self.Cbonds.perturb_bonds(getfloat("BondSDeviation", "IK_TLCSolver", self.config))
