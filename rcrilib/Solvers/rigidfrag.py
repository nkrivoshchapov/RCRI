import numpy as np
from numpy.linalg import norm
from enum import Enum

from rcrilib.Helpers import GroupItem, IK_ParameterSet, IK_Parameter, ikdof, getbool, IK_Math, createLogger
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_RigidFrag")

deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class rfcase(Enum):
    VOID = 1  # 0 atoms
    ATOM = 2  # 1 atom
    BOND = 3  # 2 atoms 1 bond
    GENERAL = 4  # N atoms, N-1 bond
    UNKNOWN = 0

class IK_RigidFrag:
    def __init__(self, mypair, config):
        self.iatoms = []  # internal
        self.ibonds = []
        self.config = config

        if mypair[1] - 1 == mypair[0] + 1:  # side
            self.case = rfcase.ATOM
            self.satoms = [mypair[0] + 1]
            self.sbonds = []
        elif mypair[1] == mypair[0] + 3:
            self.case = rfcase.BOND
            self.satoms = [mypair[0] + 1, mypair[1] - 1]
            self.sbonds = [GroupItem.bondbetween(mypair[0] + 1, mypair[0] + 2)]
        else:
            self.case = rfcase.GENERAL
            self.sbonds = [GroupItem.bondbetween(mypair[0] + 1, mypair[0] + 2),
                           GroupItem.bondbetween(mypair[1] - 1, mypair[1] - 2)]
            self.satoms = [mypair[0] + 1, mypair[1] - 1]
            curatom = mypair[0] + 2
            prevatom = 0
            while not curatom == mypair[1] - 1:
                self.iatoms.append(curatom)
                if not curatom == mypair[0] + 2:
                    self.ibonds.append(GroupItem.bondbetween(curatom, curatom - 1))
                curatom += 1

        self.ebonds = [GroupItem.bondbetween(mypair[0], mypair[0] + 1),
                       GroupItem.bondbetween(mypair[1] - 1, mypair[1])]  # ends
        self.eatoms = [mypair[0], mypair[1]]
        logger.debug("[RigidFrag] My pair = "+repr(mypair))
        logger.debug("[RigidFrag] My ebonds = "+repr(self.ebonds))
        self.initPS()
        logger.debug("Representation in RF: \n" + repr(self.PS))
        self.initGeom()

    def __repr__(self):
        line = ""
        line += "Rigid fragment. Internal atoms:\n"
        for item in self.iatoms:
            line += repr(item) + "\n"
        line += "Side atoms:\n"
        for item in self.satoms:
            line += repr(item) + "\n"
        line += "End atoms:\n"
        for item in self.eatoms:
            line += repr(item) + "\n"
        return line

    def initPS(self):
        self.PS = IK_ParameterSet()
        if self.case == rfcase.BOND:
            self.setparam(self.sbonds[0])
        elif self.case == rfcase.GENERAL:
            self.setparam(self.sbonds[0])
            for bond in self.ibonds:
                self.setparam(bond)
            self.setparam(self.sbonds[1])

    def getPS(self):
        return self.PS.getPS()

    def setparam(self, bond):
        if bond['free']:
            self.PS += IK_Parameter(ikdof.CONTINUOUS, ikdof.FREE, bond.bond_getsideatoms(), tc.RIGIDFRAG, bond=bond)
        elif not bond['free'] and bond['dep']:
            self.PS += IK_Parameter(ikdof.CONTINUOUS, ikdof.DEPENDENT, bond.bond_getsideatoms(), tc.RIGIDFRAG, bond=bond)
        elif not bond['free'] and not bond['dep']:
            self.PS += IK_Parameter(ikdof.CONTINUOUS, ikdof.FIXED, bond.bond_getsideatoms(), tc.RIGIDFRAG, bond=bond)

    def initGeom(self):
        # BOND LENGTHS
        self.elen = [self.ebonds[0].bond_getlength(), self.ebonds[1].bond_getlength()]
        self.ilen = []
        if self.case == rfcase.ATOM:
            self.slen = []
        elif self.case == rfcase.BOND:
            self.slen = [self.sbonds[0].bond_getlength()]
        elif self.case == rfcase.GENERAL:
            self.slen = [self.sbonds[0].bond_getlength(), self.sbonds[1].bond_getlength()]
            for bond in self.ibonds:
                self.ilen.append(bond.bond_getlength())

        # VAL ANGLES
        self.ival = []
        if self.case == rfcase.ATOM:
            self.sval = [self.satoms[0].atom_getangle()]
        elif self.case == rfcase.BOND:
            self.sval = [self.satoms[0].atom_getangle(), self.satoms[1].atom_getangle()]
        elif self.case == rfcase.GENERAL:
            self.sval = [self.satoms[0].atom_getangle(), self.satoms[1].atom_getangle()]
            for at in self.iatoms:
                self.ival.append(at.atom_getangle())

        # INITIALIZE FIXED PARAMETERS
        for par in self.PS:
            if par.isFixed():
                par.setValue(par.bond.bond_gettorsion())

    def applyPS(self):
        if self.case == rfcase.ATOM:
            self.leng = 0
            self.vang1 = 90
            self.vang2 = 90
            self.tang = self.satoms[0].atom_getangle() * rad2deg
        elif self.case == rfcase.BOND:
            self.leng = self.sbonds[0].bond_getlength()
            self.vang1 = self.satoms[0].atom_getangle() * rad2deg
            self.vang2 = self.satoms[1].atom_getangle() * rad2deg
            self.tang = self.PS[0].value * rad2deg
        elif self.case == rfcase.GENERAL:
            RTmat = IK_Math.Vmat(self.satoms[0].atom_getangle())@\
                    IK_Math.Dmat(self.sbonds[0].bond_getlength())@\
                    IK_Math.Tmat(self.PS[0].value)
            logger.debug("Setting torsions:")
            for item in self.PS:
                logger.debug(repr(item.value)+"("+repr(item.value*rad2deg)+")")
            logger.debug("My bonds side=%s ends=%s inter=%s"%(repr(self.sbonds),repr(self.ebonds),repr(self.ibonds)))
            self.eatoms[0]['fragframe'] = np.array([-self.ebonds[0].bond_getlength(),0.0,0.0])
            self.satoms[0]['fragframe'] = np.array([0.0, 0.0, 0.0])

            i=0
            firstside_xyz = np.array([0,0,0])
            firstside_dir = np.array([1.0,0,0])
            logger.debug("Size of PS = "+str(len(self.PS.params)))
            for atom in self.iatoms:
                atom['frame_coord'] = IK_Math.xyz_from_frame(RTmat)
                self.iatoms[i]['fragframe'] = np.array([atom['frame_coord'][0],
                                                        atom['frame_coord'][1],
                                                        atom['frame_coord'][2]])
                RTmat = RTmat@IK_Math.Vmat(atom.atom_getangle()) @ \
                              IK_Math.Dmat(atom['bonds'][1].bond_getlength()) @ \
                              IK_Math.Tmat(self.PS[i+1].value)
                logger.debug("PS[%d] was used"%(i))
                i += 1
            self.satoms[1]['fragframe'] = np.array([IK_Math.xyz_from_frame(RTmat)[0],
                                                    IK_Math.xyz_from_frame(RTmat)[1],
                                                    IK_Math.xyz_from_frame(RTmat)[2]])
            RTmat = RTmat@IK_Math.Vmat(self.satoms[1].atom_getangle())
            lastside_xyz = IK_Math.xyz_from_frame(RTmat)
            lastside_dir = IK_Math.dir_from_frame(RTmat)

            RTmat = RTmat@IK_Math.Dmat(self.ebonds[1].bond_getlength())
            self.eatoms[1]['fragframe'] = np.array([IK_Math.xyz_from_frame(RTmat)[0],
                                                    IK_Math.xyz_from_frame(RTmat)[1],
                                                    IK_Math.xyz_from_frame(RTmat)[2]])
            self.leng = norm(firstside_xyz - lastside_xyz)
            self.vang1 = np.arccos(
                np.dot(firstside_xyz - lastside_xyz, firstside_dir) / norm(firstside_xyz - lastside_xyz) / norm(
                    firstside_dir)) * rad2deg
            self.vang2 = np.arccos(
                np.dot(firstside_xyz - lastside_xyz, lastside_dir) / norm(firstside_xyz - lastside_xyz) / norm(
                    lastside_dir)) * rad2deg
            self.tang = IK_Math.gettorsion([np.array(firstside_xyz - firstside_dir), firstside_xyz, lastside_xyz,
                                            np.array(lastside_xyz + lastside_dir)]) * rad2deg

            av = firstside_dir
            bv = firstside_xyz - lastside_xyz
            cv = IK_Math.gs_rand(av, bv)

            for atom in self.iatoms:
                atom['frame_coord'] = np.array([np.dot(atom['frame_coord'], av),
                                                np.dot(atom['frame_coord'], bv),
                                                np.dot(atom['frame_coord'], cv)])
            if getbool("DoValidation", "IK_RigidFrag", self.config):
                logger.info("Running geometry check...")
                passed = True
                for atom in self.satoms+self.iatoms:
                    logger.debug("checking atom " + repr(atom))
                    if not abs(atom.atom_getangle(forcecalc=True,attr="fragframe") -
                               atom.atom_getangle(forcecalc=False,attr="fragframe")) < 0.001:
                        passed = False
                        logger.error("VAngle %s constraint is not satisfied - %f instead of %f" % (repr(atom),
                                               atom.atom_getangle(forcecalc=True,attr="fragframe"),
                                               atom.atom_getangle(forcecalc=False,attr="fragframe")))

                for bond in self.ebonds+self.ibonds+self.sbonds:
                    if not abs(bond.bond_getlength(forcecalc=True,attr="fragframe") -
                               bond.bond_getlength(forcecalc=False,attr="fragframe")) < 0.001:
                        passed = False
                        logger.error("Bondlength %s constraint is not satisfied - %f instead of %f" % (repr(bond),
                                                bond.bond_getlength(forcecalc=True,attr="fragframe"),
                                                bond.bond_getlength(forcecalc=False,attr="fragframe")))

                for param in self.PS:
                    if not abs(param.bond.bond_gettorsion(attr="fragframe") - param.value) < 0.001:
                        passed = False
                        logger.error("Torsion on %s constraint is not satisfied - %f instead of %f" % (repr(param.bond),
                                                param.bond.bond_gettorsion(attr="fragframe"),
                                                param.value))
                if not passed:
                    raise Exception("Fragprep. Not a solution")
                else:
                    logger.info("Fragment test passed")

    def get_tangle(self):
        return self.tang

    def get_vangles(self):
        return self.vang1, self.vang2

    def get_length(self):
        return self.leng

    def printgeom(self):
        eat = self.eatoms[1]
        at = self.eatoms[0]
        while True:
            print("C %14.6f %14.6f %14.6f" % (at['xyz'][0], at['xyz'][1], at['xyz'][2]))
            at += 1
            if at == eat:
                break

    def setgeom(self):
        if not self.case == rfcase.GENERAL:
            return True
        end0vec = self.eatoms[0]['xyz']
        side0vec = self.satoms[0]['xyz']
        side1vec = self.satoms[1]['xyz']

        av = side0vec - end0vec
        bv = side0vec - side1vec
        cv = IK_Math.gs_rand(av, bv)

        for atom in self.iatoms:
            atom['xyz'] = self.satoms[0]['xyz'] + atom['frame_coord'][0] * av + \
                          atom['frame_coord'][1] * bv + \
                          atom['frame_coord'][2] * cv
        if getbool("DoValidation", "IK_RigidFrag", self.config):
            logger.info("Running geometry check...")
            passed = True
            for atom in self.satoms + self.iatoms:
                logger.debug("checking atom " + repr(atom))
                if not abs(atom.atom_getangle(forcecalc=True) - atom.atom_getangle(forcecalc=False)) < 0.001:
                    passed = False
                    logger.error("VAngle %s constraint is not satisfied - %f instead of %f" % (repr(atom),
                                                        atom.atom_getangle(forcecalc=True),
                                                        atom.atom_getangle(forcecalc=False)))

            for bond in self.ebonds:
                if not abs(bond.bond_getlength(forcecalc=True) - bond.bond_getlength(forcecalc=False)) < 0.001:
                    passed = False
                    logger.error("End Bondlength %s constraint is not satisfied - %f instead of %f" % (repr(bond),
                                                        bond.bond_getlength(forcecalc=True),
                                                        bond.bond_getlength(forcecalc=False)))
            for bond in self.sbonds:
                if not abs(bond.bond_getlength(forcecalc=True) - bond.bond_getlength(forcecalc=False)) < 0.001:
                    passed = False
                    logger.error("Side Bondlength %s constraint is not satisfied - %f instead of %f" % (repr(bond),
                                                        bond.bond_getlength(forcecalc=True),
                                                        bond.bond_getlength(forcecalc=False)))
            for bond in self.ibonds:
                if not abs(bond.bond_getlength(forcecalc=True) - bond.bond_getlength(forcecalc=False)) < 0.001:
                    passed = False
                    logger.error("Internal Bondlength %s constraint is not satisfied - %f instead of %f" % (repr(bond),
                                                        bond.bond_getlength(forcecalc=True),
                                                        bond.bond_getlength(forcecalc=False)))

            for param in self.PS:
                if not abs(param.bond.bond_gettorsion()-param.value) < 0.001:
                    passed = False
                    logger.error("Torsion on %s constraint is not satisfied - %f instead of %f" % (repr(param.bond),
                                                        param.bond.bond_gettorsion(),
                                                        param.value))

            if not passed:
                raise Exception("[RIGIDFRAG] Not a solution")
            else:
                logger.info("Fragment test passed")
