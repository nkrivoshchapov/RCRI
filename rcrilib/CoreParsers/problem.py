import networkx as nx
from copy import deepcopy
from enum import Enum

from rcrilib.Solvers import IK_FourFiveCycleSolver, IK_FlappingSolver, IK_IdentitySolver, IK_Solver, IK_TLCsolver
from rcrilib.Helpers import IK_GeomValidator, createLogger, IK_ParameterSet, getbool, cause, IK_SolutionCounter, \
    sharedValue

logger = createLogger("IK_Problem")

class ikalg(Enum):
    TLC = 1
    IDENTITY = 2
    FOURFIVE = 3
    FRAGFLAP = 4

class IK_Problem:
    def __init__(self, mygraph, config, consbonds=[]):
        self.G = mygraph
        self.config = config
        self.consbonds = consbonds

        self.ndof = self.getdof()
        logger.debug("My NDOF = %d" % (self.ndof))
        logger.debug("My consbonds = " + repr(consbonds))
        self.method = self.get_solver()

    def get_solver(self):
        if IK_TLCsolver.canapply(self.G, self.ndof, self.consbonds):
            return ikalg.TLC
        elif IK_FourFiveCycleSolver.canapply(self.G, self.ndof, self.consbonds, self.config):
            IK_Solver.add_constraints(self.G)
            return ikalg.FOURFIVE
        elif IK_FlappingSolver.canapply(self.G, self.ndof, self.consbonds, self.config):
            IK_Solver.add_constraints(self.G)
            return ikalg.FRAGFLAP
        else:
            IK_Solver.add_constraints(self.G)
            return ikalg.IDENTITY

    def recheck_method(self):
        oldmethod = self.method
        self.method = self.get_solver()
        if oldmethod == ikalg.TLC and oldmethod != self.method:
            IK_Solver.add_constraints(self.G)
        return oldmethod == self.method

    def construct_solver(self):
        if self.method == ikalg.TLC:
            self.solver = IK_TLCsolver(self.G, self.ndof, self.consbonds, self.config)
        elif self.method == ikalg.FOURFIVE:
            self.solver = IK_FourFiveCycleSolver(self.G, self.ndof, self.consbonds, self.config)
        elif self.method == ikalg.FRAGFLAP:
            self.solver = IK_FlappingSolver(self.G, self.ndof, self.consbonds, self.config)
        elif self.method == ikalg.IDENTITY:
            self.solver = IK_IdentitySolver(self.G, self.ndof, self.consbonds, self.config)
        self.PS = self.solver.getPS()

        self.scounter = IK_SolutionCounter(self.PS)

        self.connectDepParams()
        if getbool("DoValidation", "IK_Problem", self.config):
            self.validator = IK_GeomValidator(self.G, self.PS)

    def getPS(self, errorcause=cause.unknown):
        return self.PS.getPS(excludeFixed=True, errorcause=errorcause)

    def getFullPS(self):
        return self.PS.getPS()

    def getdof(self):
        ncycles = nx.number_of_edges(self.G) - nx.number_of_nodes(self.G) + 1
        nnodes = len(self.G.nodes())
        solut_constr = 0
        for bond in self.consbonds:
            if len(bond) > 2:
                solut_constr += 1
            else:
                solut_constr += len(bond) - 1
        return nnodes - 3 - 3 * ncycles - solut_constr

    def connectDepParams(self):
        # TODO Create data structures for this stuff?
        # TODO Double checking the values of several LB for the same parameter (debugging);
        # TODO Remove LBs for the same parameter (release)
        for lb in self.consbonds:
            for param in self.PS:
                if (param.isDependent() or param.isFixed()) and len(param.atoms) == 2:
                    bond = [param.atoms[0], param.atoms[1]]
                    if lb.hasBond(bond):
                        if not lb.isFake():
                            lb.setParam(param)
                        param['linkingbond'] = lb
                        lb.checkOrientation()
                        logger.debug("Bond of dependent parameter = " + repr(bond))
                        break

    def applyPS(self):
        if getbool("SmartAPS", "IK_Problem", self.config) and not self.PS.isChanged():
            logger.info("Skipping solver run. LS = " + repr(self.PS.cause))
            if self.PS.cause == cause.success:
                return IK_ParameterSet()
            else:
                return self.PS.getPS(excludeFixed=True, errorcause=self.PS.cause)

        redoPS = self.solver.applyPS()
        if redoPS.success:
            self.loadGraphCoord()
            if self.solver == ikalg.TLC:
                redoPS = self.solver.checkoverlap()
        logger.debug("redoPS = " + repr(redoPS))
        if redoPS.success:
            self.PS.cause = cause.success
        else:
            self.PS.cause = redoPS.cause
        return redoPS.getPS(errorcause=redoPS.cause)

    def loadGraphCoord(self):
        self.solver.updateGraphxyz()
        if getbool("DoValidation", "IK_Problem", self.config):
            good = self.validator.validate(self.G)
            logger.info("Validator returned "+repr(good))
            if not good:
                raise Exception("[Problem] Failed to obtain IK solution")

    def updatexyz(self, G):
        self.solver.updatexyz(G)

    def perturb_geometry(self):
        if self.method == ikalg.TLC:
            self.solver.perturb_geometry()

    def share_vangle(self, atom):
        if "shared_vangle" not in self.G.nodes[atom]:
            self.G.nodes[atom]['own_vangle'] = True
            self.G.nodes[atom]['shared_vangle'] = sharedValue("Shared vangle on atom %s (nb = %s)" % (atom,
                                                                              repr(list(self.G.neighbors(atom)))))
        return self.G.nodes[atom]['shared_vangle']

    def assign_shared_vangle(self, atom, svalue):
        if "shared_vangle" in self.G.nodes[atom]:
            myvalue = self.G.nodes[atom]['shared_vangle']
            myvalue.value = svalue
        else:
            self.G.nodes[atom]['shared_vangle'] = svalue
        self.G.nodes[atom]['own_vangle'] = False

    def share_length(self, bond):
        if "shared_length" not in self.G[bond[0]][bond[1]]:
            self.G[bond[0]][bond[1]]['own_length'] = True
            self.G[bond[0]][bond[1]]['shared_length'] = sharedValue("Shared bond length (%d,%d)" % (bond[0], bond[1]))
        return self.G[bond[0]][bond[1]]['shared_length']

    def assign_shared_length(self, bond, svalue):
        if "shared_length" in self.G[bond[0]][bond[1]]:
            myvalue = self.G[bond[0]][bond[1]]['shared_length']
            myvalue.value = svalue
        else:
            self.G[bond[0]][bond[1]]['shared_length'] = svalue
        self.G[bond[0]][bond[1]]['own_length'] = False
