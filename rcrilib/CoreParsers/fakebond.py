import networkx as nx
from copy import deepcopy

from .linkingbond import IK_LinkingBond
from rcrilib.Helpers import Segment, itemtype, createLogger, IK_ParameterSet, IK_Parameter, ikdof, IK_Math, getbool, \
                            IK_GeomValidator
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_FakeBond")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class IK_FakeBond:
    def __init__(self, nodegraph, nbgraph, atoms, ind_cycle_idx, dep_cycle_idx, config):
        self.cycleG = nodegraph  # saving a pointer to graph in IK_CyclicPart.IG
        self.atom_num = len(atoms) + 2
        self.bond_num = len(atoms) + 1
        self.atoms = Segment(self.atom_num, itype=itemtype.ATOM)
        self.bonds = Segment(self.bond_num, itype=itemtype.BOND)
        self.config = config

        tempgraph = nx.Graph()
        for edge in self.cycleG.edges(data = True):
            if edge[0] in atoms or edge[1] in atoms:
                tempgraph.add_edge(edge[0], edge[1])
                for attr in edge[2]:
                    tempgraph[edge[0]][edge[1]][attr] = edge[2][attr]
        for node in self.cycleG.nodes(data = True):
            if tempgraph.has_node(node[0]):
                for attr in node[1]:
                    tempgraph.nodes[node[0]][attr] = node[1][attr]

        nextnode = []
        for node in tempgraph:
            if len(list(tempgraph.neighbors(node))) == 1:
                nextnode.append(node)
                break
        prevnode = -1
        i = 0
        while len(nextnode) > 0:
            if len(nextnode) > 1:
                raise Exception("Error")
            mynode = nextnode[0]
            self.atoms[i]['G_idx'] = mynode
            nextnode = list(tempgraph.neighbors(mynode))
            if prevnode != -1:
                nextnode.remove(prevnode)
            prevnode = mynode
            i += 1

        self.atoms.initGidx()
        self.atoms.set_bond_seg(self.bonds)
        self.bonds.set_atom_seg(self.atoms)
        for i in range(self.bond_num):
            self.atoms[i]['bond'] = self.bonds[i]
            self.bonds[i]['atoms'] = [self.atoms[i], self.atoms[i+1]]

        for node in tempgraph.nodes(data=True):
            for attr in node[1]:
                self.atoms.getbyGidx(node[0])[attr] = node[1][attr]
        for edge in tempgraph.edges(data = True):
            for attr in edge[2]:
                self.bonds.getbyGidx([edge[0],edge[1]])[attr] = edge[2][attr]

        for bond in self.bonds.getinner():
            depbond = [bond['atoms'][0]['G_idx'], bond['atoms'][1]['G_idx']]
            node_sideat = [list(nodegraph.neighbors(depbond[0])), list(nodegraph.neighbors(depbond[1]))]
            nb_sideat = [list(nbgraph.neighbors(depbond[0])), list(nbgraph.neighbors(depbond[1]))]
            for item in [node_sideat[0], nb_sideat[0]]:
                item.remove(depbond[1])
            for item in [node_sideat[1], nb_sideat[1]]:
                item.remove(depbond[0])
            logger.info("Linking bond: %s ; Side1: %s ; Side2 : %s" % (repr(depbond),
                                                                repr([node_sideat[0][0], node_sideat[1][0]]),
                                                                repr([nb_sideat[0][0], nb_sideat[1][0]])))
            bond['lb_inst'] = IK_LinkingBond(depbond, [node_sideat[0][0], node_sideat[1][0]],
                                             [nb_sideat[0][0], nb_sideat[1][0]],
                                             ind_cycle_idx, dep_cycle_idx)

        logger.info("My nodes: " + repr(list(self.cycleG.nodes())))
        logger.info("My edges: " + repr(list(self.cycleG.edges())))

        self.PS = IK_ParameterSet()
        for bond in self.bonds.getinner():
             newpar = IK_Parameter(ikdof.CONTINUOUS, ikdof.DEPENDENT, bond.bond_getsideatoms(), tc.FAKEBOND,
                                    True, bond=bond)
             bond['lb_inst'].setParam(newpar)
             bond['lb_inst'].checkOrientation()
             newpar.makeshared("Value of dihedral (%d,%d,%d,%d) on >2 atoms linkage" % (bond.bond_getsideatoms()[0],
                                                                                        bond['atoms'][0]['G_idx'],
                                                                                        bond['atoms'][1]['G_idx'],
                                                                                        bond.bond_getsideatoms()[1]))
             self.cycleG[bond['atoms'][0]['G_idx']] \
                        [bond['atoms'][1]['G_idx']]['shared_dihedral'] = newpar.value # Get reference, not value itself
             self.PS += newpar

        self.ind_idx = ind_cycle_idx
        self.dep_idx = dep_cycle_idx
        self.eatoms = self.atoms.getends()
        self.satoms = self.atoms.getsides()
        self.iatoms = self.atoms.getinternal()
        self.ebonds = self.bonds.getends()
        self.sbonds = self.bonds.getsides()
        self.ibonds = self.bonds.getinternal()
        self.bond = [self.satoms[0]['G_idx'], self.satoms[1]['G_idx']]
        self.inner_idx = []
        for atom in self.atoms.getinner():
            self.inner_idx.append(atom['G_idx'])

    def readDependence(self, G):
        for bond in self.bonds.getinner():
            bond['lb_inst'].readDependence(G)

    def resolveDependence(self, ind_G):
        for bond in self.bonds.getinner():
            bond['lb_inst'].resolveDependence(ind_G)

    def __eq__(self, other):
        if isinstance(other, list):
            return other[0] in self.inner_idx and other[1] in self.inner_idx
        else:
            raise Exception("Invalid comparison")

    def hasBond(self, atompair):
        if isinstance(atompair, list):
            return atompair[0] in self.inner_idx and atompair[1] in self.inner_idx
        else:
            raise Exception("Invalid comparison")

    def checkOrientation(self):
        pass # TODO Should get rid of this function at all

    def __len__(self):
        return self.atom_num

    def getPS(self):
        pass

    def isFake(self):
        return True

    def connect_cycles(self, G_dep, G_ind, G_glob):
        for atom in self.atoms:
            atom['xyz'] = deepcopy(self.cycleG.nodes[atom['G_idx']]['xyz'])

        self.bonds[1]['lb_inst'].connect_cycles(self.cycleG, G_ind, G_glob)
