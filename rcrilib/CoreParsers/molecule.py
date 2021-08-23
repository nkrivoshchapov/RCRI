import networkx as nx
from copy import deepcopy, copy
import numpy as np
from numpy.linalg import norm, inv
import pickle, random

from rcrilib.Helpers import IK_ParameterSet, IK_Parameter, ikdof, getfloat, getbool, cause, IK_GeomValidator, \
    geomcheckTree_graph, IK_Math, createLogger, IK_CompoundCounter
from rcrilib.Helpers import target_class as tc
from .cyclicpart import IK_CyclicPart

logger = createLogger("IK_Molecule")

class IK_Molecule:
    def __init__(self, sdffile, config):
        self.config = config
        self.atoms_sym = []
        self.bonds_num = []
        self.bonds_type = []
        self.readsdf(sdffile)
        self.discr_cp = -1

    def readsdf(self, file):
        lines = open(file, "r").readlines()
        natoms = int(lines[3][0:3])
        nbonds = int(lines[3][3:6])
        logger.debug("NAtoms = %d   NBonds = %d" % (natoms, nbonds))
        self.atoms_xyz = np.zeros((natoms, 3))
        for i in range(4, 4 + natoms):
            parts = lines[i].replace("\n", "").split()
            self.atoms_xyz[i-4][0] = float(parts[0])
            self.atoms_xyz[i-4][1] = float(parts[1])
            self.atoms_xyz[i-4][2] = float(parts[2])
            self.atoms_sym.append(parts[3])
        for i in range(4 + natoms, 4 + natoms + nbonds):
            parts = lines[i].replace("\n", "").split()
            at1 = int(lines[i][0:3])
            at2 = int(lines[i][3:7])
            bondtype = int(lines[i][7:10])
            self.bonds_num.append([at1 - 1, at2 - 1])
            self.bonds_type.append(bondtype)

    def buildmolgraph(self):
        self.molgr = nx.Graph()
        self.molgr.add_edges_from(self.bonds_num)
        for i in range(len(self.atoms_xyz)):
            self.molgr.nodes()[i]['xyz'] = np.array(self.atoms_xyz[i])

    def hasaneighbour(self, atom, otheratom):
        ans = False
        for bond in self.bonds_num:
            if atom in bond and otheratom not in bond:
                ans = True
                break
        return ans

    def getaneighbour(self, atom, otheratom):
        good_ones = set(self.molgr.neighbors(atom))
        good_ones.discard(otheratom)
        for cp in self.cycparts:
            if atom in cp.G.nodes():
                return good_ones.intersection(set(cp.G.nodes())).pop()
        else:
            return good_ones.pop()

    def getbondidx(self, mybond):
        for i in range(len(self.bonds_num)):
            if mybond[0] in self.bonds_num[i] and mybond[1] in self.bonds_num[i]:
                return i

    def getfreebonds(self):
        self.buildmolgraph()
        freebonds = list(nx.bridges(self.molgr))
        for bond in reversed(freebonds):
            if not (self.hasaneighbour(bond[0], bond[1]) and self.hasaneighbour(bond[1], bond[0])):
                freebonds.remove(bond)
        for bond in reversed(freebonds):
            if self.bonds_type[self.getbondidx(bond)] > 1:
                freebonds.remove(bond)
        self.freebonds = deepcopy(freebonds)

    def prepare_mc(self):
        self.getfreebonds()
        for bond in self.freebonds:
            logger.debug(repr([bond[0] + 1, bond[1] + 1]))

    def getrotbonds(self):
        self.buildmolgraph()
        rotbonds = deepcopy(self.bonds_num)
        for bond in reversed(rotbonds):
            if not (self.hasaneighbour(bond[0], bond[1]) and self.hasaneighbour(bond[1], bond[0])):
                rotbonds.remove(bond)
        for bond in reversed(rotbonds):
            if self.bonds_type[self.getbondidx(bond)] > 1:
                rotbonds.remove(bond)
        self.rotbonds = deepcopy(rotbonds)

    def gencyclparts(self):
        logger.info("Largest cycle has %d atoms" % max([len(c) for c in nx.minimum_cycle_basis(self.molgr)]))
        bridges = [list(c) for c in list(nx.bridges(self.molgr))]
        nonbridges = []
        for bond in self.bonds_num:
            if bond not in bridges and [bond[1], bond[0]] not in bridges:
                nonbridges.append(copy(bond))
        cyclgraph = nx.Graph()
        cyclgraph.add_edges_from(nonbridges)
        subgraphs = [cyclgraph.subgraph(c) for c in nx.connected_components(cyclgraph)]
        self.cycparts = []
        logger.debug("%d cyclic parts has been found in the molecule" % (len(subgraphs)))

        for gr in subgraphs:
            logger.debug("Processing cyclic part %d/%d" % (subgraphs.index(gr) + 1, len(subgraphs)))
            self.cycparts.append(IK_CyclicPart(gr, self.molgr, self.atoms_xyz,
                                               self.bonds_num, self.bonds_type, self.config))

        self.PS = IK_ParameterSet()
        self.FullPS = IK_ParameterSet()
        for part in self.cycparts:
            self.PS += part.getPS()
            self.FullPS += part.getFullPS()

        self.ccounter = IK_CompoundCounter()
        for part in self.cycparts:
            self.ccounter += part.ccounter
        for i, item in enumerate(self.ccounter.counters):
            item.name = str(i)
        self.noncyclicDOF()
        self.buildfraggraph()
        self.builddirectedgraph()
        if getbool("DoValidation", "IK_Molecule", self.config):
            self.validator = IK_GeomValidator(self.molgr, self.FullPS)
        logger.debug("Full Mol PS: " + repr(self.PS))

    def buildfraggraph(self):
        # TODO FIX BUG WHEN NO NON-CYCLIC DOFs
        fragments = deepcopy(self.molgr)
        for par in self.TorPS:
            fragments.remove_edge(par.atoms[0], par.atoms[1])
        frag_comp = [set(fragments.subgraph(c).nodes()) for c in
                   nx.connected_components(fragments)]
        fg_params = []
        fg_atoms = []
        for i, comp in enumerate(frag_comp):
            newitem = set()
            for j in comp:
                newitem.add(j)
            fg_atoms.append(newitem)

            newps = IK_ParameterSet()
            for par in self.TorPS:
                if par.atoms[0] in newitem or par.atoms[1] in newitem:
                    newps += par
                    par.frags.append(i)
            if len(newps.params) > 0:
                fg_params.append(newps)

        fg_edges = []
        for par in self.TorPS:
            fg_edges.append(deepcopy(par.frags))

        self.FG = nx.Graph()
        if len(frag_comp) == 1:
            self.FG.add_node(0)
            self.FG.nodes[0]['atoms'] = fg_atoms[0]
        else:
            self.FG.add_edges_from(fg_edges)
            for par in self.TorPS:
                self.FG[par.frags[0]][par.frags[1]]['param'] = IK_ParameterSet(p = par)

            for i in range(len(fg_atoms)):
                self.FG.nodes[i]['params'] = fg_params[i]
                self.FG.nodes[i]['atoms'] = fg_atoms[i]

        for cycpart in self.cycparts:
            atom = list(cycpart.G.nodes)[0]
            for i in range(len(self.FG.nodes)):
                if atom in self.FG.nodes[i]['atoms']:
                    self.FG.nodes[i]['cycpart'] = cycpart
                    break

        for i in range(len(self.FG.nodes)):
            if 'cycpart' not in self.FG.nodes[i]:
                self.FG.nodes[i]['cycpart'] = False

    def builddirectedgraph(self):
        self.headnode = 0
        endnodes = []
        for node in list(self.FG):
            if len(list(self.FG.neighbors(node))) == 1:
                endnodes.append(node)
        if self.headnode in endnodes:
            endnodes.remove(self.headnode)

        self.DFG = nx.DiGraph(directed=True)
        self.DFG.add_nodes_from(self.FG.nodes(data=True))
        logger.debug("My endnodes " + repr(endnodes))
        for enode in endnodes:
            solvepath = nx.shortest_path(self.FG, source=enode, target=self.headnode)
            for i in range(len(solvepath) - 1):
                self.DFG.add_edge(solvepath[i], solvepath[i + 1])
                self.DFG[solvepath[i]][solvepath[i + 1]]['param'] = self.FG[solvepath[i + 1]][solvepath[i]]['param']

        self.localcoordinates = []
        for i in range(len(self.atoms_xyz)):
            self.localcoordinates.append(np.array(self.atoms_xyz[i]))
        for i in range(len(self.DFG.nodes)):
            if not self.FG.nodes[i]['cycpart']:
                if len(list(self.DFG.neighbors(i))) == 1:
                    othernode = list(self.DFG.neighbors(i))[0]
                    bond = self.DFG[i][othernode]['param'][0].atoms
                    sides = self.DFG[i][othernode]['param'][0].sides
                    if bond[1] in self.DFG.nodes[i]['atoms']:
                        at0 = np.array(self.atoms_xyz[bond[0]])
                        at1 = np.array(self.atoms_xyz[bond[1]])
                        at2 = np.array(self.atoms_xyz[sides[1]])
                    elif bond[0] in self.DFG.nodes[i]['atoms']:
                        at0 = np.array(self.atoms_xyz[bond[1]])
                        at1 = np.array(self.atoms_xyz[bond[0]])
                        at2 = np.array(self.atoms_xyz[sides[0]])
                    xv = at1 - at0
                    yv = at2 - at1
                    zv = IK_Math.gs_rand(xv, yv)
                elif len(list(self.DFG.neighbors(i))) == 0:
                    at1 = np.zeros(3)
                    xv = np.array([1, 0, 0])
                    yv = np.array([0, 1, 0])
                    zv = np.array([0, 0, 1])
                else:
                    raise Exception("self.DFG IS NOT A TREE!!!")
                for atom in self.DFG.nodes[i]['atoms']:
                    self.localcoordinates[atom] -= at1
                    self.localcoordinates[atom] = np.array([self.localcoordinates[atom] @ xv,
                                                            self.localcoordinates[atom] @ yv,
                                                            self.localcoordinates[atom] @ zv,
                                                            1])
                newframe = np.zeros((4, 4))
                newframe[:3, 0] = xv
                newframe[:3, 1] = yv
                newframe[:3, 2] = zv
                newframe[:3, 3] = at1
                newframe[3, 3] = 1
                self.DFG.nodes[i]['frame'] = newframe

        """
        # CREATE FRAGMENT'S XYZs
        for i in range(len(self.DFG.nodes)):
            myfile = open("frag_%d.xyz" % i, "w")
            for atom in self.DFG.nodes[i]['atoms']:
                myfile.write("%2s %14.6f %14.6f %14.6f\n" % (self.atoms_sym[atom],
                                                             self.localcoordinates[atom][0],
                                                             self.localcoordinates[atom][1],
                                                             self.localcoordinates[atom][2]))
            myfile.close()
        """

        for edge in list(self.DFG.edges()):
            if not self.FG.nodes[edge[1]]['cycpart']:
                startframe = self.DFG.nodes[edge[1]]["frame"]
                endbond = self.DFG[edge[0]][edge[1]]['param'][0].atoms
                endsides = self.DFG[edge[0]][edge[1]]['param'][0].sides
                if endbond[1] in self.DFG.nodes[edge[1]]['atoms']:
                    at0 = np.array(self.atoms_xyz[endbond[0]])
                    at1 = np.array(self.atoms_xyz[endbond[1]])
                    at2 = np.array(self.atoms_xyz[endsides[1]])
                elif endbond[0] in self.DFG.nodes[edge[1]]['atoms']:
                    at0 = np.array(self.atoms_xyz[endbond[1]])
                    at1 = np.array(self.atoms_xyz[endbond[0]])
                    at2 = np.array(self.atoms_xyz[endsides[0]])
                xv = at0 - at1
                yv = at2 - at1
                zv = IK_Math.gs_rand(xv, yv)
                newframe = np.zeros((4, 4))
                newframe[:3, 0] = xv
                newframe[:3, 1] = yv
                newframe[:3, 2] = zv
                newframe[:3, 3] = at0
                newframe[3, 3] = 1
                self.DFG[edge[0]][edge[1]]['transition'] = inv(startframe) @ newframe

        for i in range(len(self.FG.nodes)):
            if self.FG.nodes[i]['cycpart']:
                self.FG.nodes[i]['cycpart'].setconfig(self.molgr, self.DFG, i)

        self.DFG = self.DFG.reverse(copy=False)

    def noncyclicDOF(self):
        self.TorPS = IK_ParameterSet()
        for bridge in nx.bridges(self.molgr):
            if (self.hasaneighbour(bridge[0], bridge[1]) and self.hasaneighbour(bridge[1], bridge[0])):
                sideatoms = [self.getaneighbour(bridge[0], bridge[1]), self.getaneighbour(bridge[1], bridge[0])]
                newparam = IK_Parameter(ikdof.CONTINUOUS, ikdof.FREE, sideatoms, tc.MOLECULE, atoms=bridge)
                self.molgr[bridge[0]][bridge[1]]['par'] = newparam
                self.PS += newparam
                self.TorPS += newparam

    def __repr__(self):
        resstring = str(len(self.atoms_xyz)) + "\n\n"
        for i in range(len(self.atoms_xyz)):
            resstring += "%2s %14.6f %14.6f %14.6f\n" % (self.atoms_sym[i],
                                                         self.atoms_xyz[i][0],
                                                         self.atoms_xyz[i][1],
                                                         self.atoms_xyz[i][2])
        return resstring

    def getPS(self):
        return self.PS.getPS()

    def startDiscreteRun(self):
        self.discr_cp = 0
        for part in self.cycparts:
            part.startDiscreteRun()
        self.ccounter.startDiscreteRun()
        self.ccounter.backup_ddof()

    def done_all_discr(self):
        for part in self.cycparts:
            if not part.done_all_discr():
                return False
        return True

    def applyPS(self, increase_discrete = False):
        if increase_discrete:
            self.discr_cp = 0
            while self.cycparts[self.discr_cp].done_all_discr():
                self.cycparts[self.discr_cp].resetDDOF()
                self.discr_cp += 1
        if increase_discrete:
            redoPS = self.cycparts[self.discr_cp].applyPS(increase_discrete=increase_discrete)
            self.ccounter.checkPS(redoPS)
            logger.debug("[molecule] redoPS = " + repr(redoPS))
            if not redoPS.success:
                return redoPS
        for i in range(len(self.cycparts)):
            if i < self.discr_cp:
                self.cycparts[i].ccounter.restore_ddof()
            if increase_discrete and i == self.discr_cp:
                continue
            else:
                redoPS = self.cycparts[i].applyPS()
            self.ccounter.checkPS(redoPS)
            logger.debug("[molecule] redoPS = " + repr(redoPS))
            if not redoPS.success:
                return redoPS
            if i < self.discr_cp:
                self.cycparts[i].startDiscreteRun()

        for i in range(len(self.FG.nodes)):
            if self.FG.nodes[i]['cycpart']:
                self.FG.nodes[i]['cycpart'].updateattributes(self.DFG, i, self.localcoordinates)

        todo_frags = [self.headnode]
        while len(todo_frags) > 0:
            curfrag = todo_frags.pop()
            for atom in self.DFG.nodes[curfrag]['atoms']:
                self.atoms_xyz[atom] = list((self.DFG.nodes[curfrag]['frame'] @ self.localcoordinates[atom])[:3])
            dep_frags = list(self.DFG.neighbors(curfrag))
            for frag in dep_frags:
                self.DFG.nodes[frag]['frame'] = self.DFG.nodes[curfrag]['frame'] @ \
                                                self.DFG[curfrag][frag]['transition'] @ \
                                                IK_Math.Tmat(self.DFG[curfrag][frag]['param'][0].value)
            todo_frags += dep_frags

        for i in range(len(self.atoms_xyz)):
            self.molgr.nodes[i]['xyz'] = np.array(self.atoms_xyz[i])

        if getbool("DoValidation", "IK_Molecule", self.config):
            good = self.validator.validate(self.molgr)
            logger.info("Validator returned " + repr(good))
            if not good:
                logger.error("[Molecule] Failed to obtain IK solution")

        if (not getbool("RunGeomCheck", "IK_Molecule", self.config)) or \
                getbool("RunGeomCheck", "IK_Molecule", self.config) and \
                geomcheckTree_graph(self.molgr, getfloat("MinDistance", "IK_Molecule", self.config)):
            psout = IK_ParameterSet()
        else:
            if not self.PS.byTarget(tc.MOLECULE).isEmpty():
                psout = self.PS.byTarget(tc.MOLECULE, errorcause=cause.geomoverlap)
            else:
                doflog = open("doflog.csv", "a")
                doflog.write("Returned full PS\n")
                doflog.close()
                psout = self.PS.getPS(errorcause=cause.geomoverlap)
        return psout

    def perturb_geometry(self):
        for cp in self.cycparts:
            cp.perturb_geometry()

    def prepare_ik(self):
        self.getrotbonds()
        self.gencyclparts()

    def writeToXyz(self, name):
        lines = [str(len(self.atoms_xyz)), ""]
        for i in range(len(self.atoms_xyz)):
            lines.append("%2s %14.6f %14.6f %14.6f" % (self.atoms_sym[i],
                                                         self.atoms_xyz[i][0],
                                                         self.atoms_xyz[i][1],
                                                         self.atoms_xyz[i][2]))
        file = open(name, "w")
        file.write("\n".join(lines))
        file.close()

    def writeToMol(self, name):
        lines = ["", "", ""]
        lines.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (len(self.atoms_xyz), len(self.bonds_num)))
        #   -5.5250    1.6470    1.0014 C   0  0  0  0  0  0  0  0  0  0  0  0
        for i in range(len(self.atoms_xyz)):
            lines.append("%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0" % (
                                                       self.atoms_xyz[i][0],
                                                       self.atoms_xyz[i][1],
                                                       self.atoms_xyz[i][2],
                                                       self.atoms_sym[i]))

        for i in range(len(self.bonds_num)):
            lines.append("%3s%3s%3s  0" % (self.bonds_num[i][0] + 1,
                                           self.bonds_num[i][1] + 1,
                                           self.bonds_type[i]))
        lines.append("M  END\n")

        file = open(name, "w")
        file.write("\n".join(lines))
        file.close()

    def geomcheckGraph(self):
        return geomcheckTree_graph(self.molgr, getfloat("MinDistance", "IK_Molecule", self.config))

    def geomcheckXyz(self):
        for i in range(len(self.atoms_xyz)):
            for j in range(i + 1, len(self.atoms_xyz)):
                if [i, j] not in self.bonds_num and [j, i] not in self.bonds_num and \
                        (self.atoms_xyz[i][0] - self.atoms_xyz[j][0]) ** 2 + \
                        (self.atoms_xyz[i][1] - self.atoms_xyz[j][1]) ** 2 + \
                        (self.atoms_xyz[i][2] - self.atoms_xyz[j][2]) ** 2 < \
                        getfloat("MinDistance", "IK_Molecule", self.config)**2:  # SQUARE!
                    logger.info("Overlap detected")
                    return False
        return True
