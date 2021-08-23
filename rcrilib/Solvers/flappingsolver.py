import numpy as np
from copy import copy, deepcopy
import networkx as nx
import random
from numpy import sin, cos
from numpy.linalg import norm, det, inv

from .solver import IK_Solver
from rcrilib.Helpers import IK_ParameterSet, IK_Parameter, ikdof, getint, createLogger, IK_Math
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_FlappingSolver")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class IK_FlappingSolver(IK_Solver):
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

    @staticmethod
    def get_other_nb(nblist, c_vert):
        if nblist[0] == c_vert:
            return nblist[1]
        elif nblist[1] == c_vert:
            return nblist[0]
        else:
            raise Exception("redG has unexpected topology")

    @staticmethod
    def get_appropriate_pairs(G, config):
        redG = nx.MultiGraph()
        for edge in G.edges():
            redG.add_edge(edge[0], edge[1], key=0, path_verticies=[])
        branching_vert = []
        for node in redG:
            if len(list(redG.neighbors(node))) > 2:
                branching_vert.append(node)
        angle_threshold = getint("AngleThreshold", "IK_FlappingSolver", config)
        suppG = nx.Graph()
        for c_vert in branching_vert:
            nb_set = set(redG.neighbors(c_vert))
            nb_set.difference_update(branching_vert)
            while len(nb_set) > 0:
                p_vert = nb_set.pop()
                if len(redG.get_edge_data(c_vert, p_vert, key=0)['path_verticies']) == 0:
                    redG.get_edge_data(c_vert, p_vert, key=0)['path_verticies'].append(c_vert)
                othernb = IK_FlappingSolver.get_other_nb(list(redG.neighbors(p_vert)), c_vert)
                if len(redG.get_edge_data(p_vert, othernb, key=0)['path_verticies']) > 0:
                    raise Exception("G reduction took an unexpected turn")

                if othernb not in branching_vert:
                    redG.add_edge(c_vert, othernb, key=0)
                    redG.get_edge_data(c_vert, othernb, key=0)['path_verticies'] = \
                        redG.get_edge_data(c_vert, p_vert, key=0)['path_verticies'] + [p_vert]
                    redG.remove_node(p_vert)
                else:
                    if not suppG.has_edge(c_vert, othernb):
                        suppG.add_edge(c_vert, othernb)
                        suppG[c_vert][othernb]['path_count'] = 1
                    else:
                        suppG[c_vert][othernb]['path_count'] += 1
                    redG.add_edge(c_vert, othernb, key=suppG[c_vert][othernb]['path_count'] - 1)
                    redG.get_edge_data(c_vert, othernb,
                                            key=suppG[c_vert][othernb]['path_count'] - 1)['path_verticies'] = \
                        redG.get_edge_data(c_vert, p_vert, key=0)['path_verticies'] + [p_vert] + [othernb]
                    redG.remove_node(p_vert)
                nb_set.update(redG.neighbors(c_vert))
                nb_set.difference_update(branching_vert)

        appropriate_pairs = []
        for edge in suppG.edges():
            for c_key in range(suppG.get_edge_data(*edge)['path_count']):
                c_path = redG.get_edge_data(*edge, key=c_key)['path_verticies']
                for i in range(len(c_path) - 2):
                    a_edge = (c_path[i], c_path[i + 1])
                    for j in range(i + 3, len(c_path) - 1):
                        b_edge = (c_path[j], c_path[j + 1])
                        logger.debug("Trying pair (%d,%d) & (%d,%d)" % (*a_edge, *b_edge))
                        ab_torsion = abs(IK_Math.gettorsion([G.nodes[a_edge[0]]['xyz'],
                                                             G.nodes[a_edge[1]]['xyz'],
                                                             G.nodes[b_edge[0]]['xyz'],
                                                             G.nodes[b_edge[1]]['xyz']]))
                        logger.debug("Torsion of pair (%d,%d) & (%d,%d) = %f deg" % (*a_edge, *b_edge,
                                                                                     ab_torsion * rad2deg))
                        if ab_torsion < angle_threshold * deg2rad:
                            logger.debug("Pair (%d,%d) & (%d,%d) passed torsion check" % (*a_edge, *b_edge))
                            middle_atoms = tuple(c_path[l] for l in range(i + 2, j))
                            appropriate_pairs.append(a_edge + middle_atoms + b_edge)
                        # TODO Use 'single' edge attribute to exclude distortion of double,triple bonds
        return appropriate_pairs

    def genconformers(self):
        self.appropriate_pairs = IK_FlappingSolver.get_appropriate_pairs(self.G, self.config)
        # TODO Implement the check for independency between self.appropriate_pairs, to remove redundant DDOFs
        # TODO Geom Check at this stage
        # TODO Scalable conformer generation OR Implement both cached mode and it-time mode
        # TODO Smaller conformer set when 2**len(ap) is too much (max K flipping steps)
        self.N = 2 ** len(self.appropriate_pairs)
        logger.info("We'll have %d flipping combinations!" % self.N)
        self.nodes = list(self.G.nodes())
        if self.config["IK_FlappingSolver"].get("GenerationMode") == "intime":
            self.conformer = np.empty(shape=(len(self.G.nodes()), 3))
            self.starting_geom = np.empty(shape=(len(self.G.nodes()), 3))
            for atom_idx in range(len(self.nodes)):
                self.starting_geom[atom_idx][:] = self.G.nodes[self.nodes[atom_idx]]['xyz'][:]
        else:
            self.conformers = np.empty(shape=(self.N, len(self.G.nodes()), 3))
            for conf_idx in range(self.N):
                for atom_idx in range(len(self.nodes)):
                    self.conformers[conf_idx][atom_idx][:] = self.G.nodes[self.nodes[atom_idx]]['xyz'][:]
            for idx_frag, c_frag in enumerate(self.appropriate_pairs):
                logger.info("Processing pair %d/%d" % (idx_frag, len(self.appropriate_pairs)))
                for conf_idx in range(self.N):
                    if conf_idx % 2**(idx_frag + 1) >= 2**idx_frag:
                        at0 = self.conformers[conf_idx][self.nodes.index(c_frag[0])]
                        at1 = self.conformers[conf_idx][self.nodes.index(c_frag[1])]
                        at2 = self.conformers[conf_idx][self.nodes.index(c_frag[len(c_frag) - 2])]
                        at3 = self.conformers[conf_idx][self.nodes.index(c_frag[len(c_frag) - 1])]

                        av = at2 - at1
                        av /= norm(av)
                        bv = at0 - at1
                        cv = np.cross(av, bv)
                        cv /= norm(cv)
                        bv = np.cross(cv, av)

                        c_frame = np.array([av, bv, cv])
                        tordiff = IK_Math.gettorsion_vec(bv, av, at3 - at2)
                        corr_matrix = np.array([[1, 0, 0],
                                                [0, cos(tordiff/2), sin(tordiff/2)],
                                                [0, -sin(tordiff/2), cos(tordiff/2)]])
                        corr_matrix = inv(c_frame) @ corr_matrix @ c_frame
                        bv = corr_matrix @ bv
                        cv = corr_matrix @ cv
                        c_frame = np.array([av, bv, cv])
                        reflection_matrix = np.array([[1, 0, 0],
                                                      [0, 1, 0],
                                                      [0, 0, -1]])
                        reflection_matrix = inv(c_frame) @ reflection_matrix @ c_frame
                        logger.debug(repr(det(reflection_matrix)))

                        for atom in c_frag[2:len(c_frag) - 2]:
                            self.conformers[conf_idx][self.nodes.index(atom)][:] = at1 + reflection_matrix @ \
                                                              (self.conformers[conf_idx][self.nodes.index(atom)] - at1)
        self.PS.getMyDiscreteParameter(tc.SOLVER).maxValue = self.N

    @staticmethod
    def canapply(G, ndof, consbonds, config):
        # Can apply if there is something to flip
        Npairs = len(IK_FlappingSolver.get_appropriate_pairs(G, config))
        if Npairs > 9 and config["IK_FlappingSolver"].get("GenerationMode") != "intime":
            logger.error("%d conformers is too much to Handel! Initialization will take long! Use \"intime\" key "
                         "or decrease AngleThreshold" % 2 ** Npairs)
        return Npairs > 0

    def updatexyz(self, G):
        logger.info("Getting conformer #%d" % self.PS.getMyDiscreteParameter(tc.SOLVER).value)
        if self.config["IK_FlappingSolver"].get("GenerationMode") == "intime":
            for atom_idx in range(len(self.nodes)):
                G.nodes[self.nodes[atom_idx]]['xyz'][:] = self.conformer[atom_idx][:]
        else:
            for atom_idx in range(len(self.nodes)):
                G.nodes[self.nodes[atom_idx]]['xyz'][:] = \
                    self.conformers[self.PS.getMyDiscreteParameter(tc.SOLVER).value][atom_idx][:]

    def updateGraphxyz(self):
        logger.info("Getting conformer #%d" % self.PS.getMyDiscreteParameter(tc.SOLVER).value)
        if self.config["IK_FlappingSolver"].get("GenerationMode") == "intime":
            for atom_idx in range(len(self.nodes)):
                self.G.nodes[self.nodes[atom_idx]]['xyz'][:] = self.conformer[atom_idx][:]
        else:
            for atom_idx in range(len(self.nodes)):
                self.G.nodes[self.nodes[atom_idx]]['xyz'][:] = \
                    self.conformers[self.PS.getMyDiscreteParameter(tc.SOLVER).value][atom_idx][:]

    def getPS(self):
        return self.PS.getPS()

    def applyPS(self):
        # TODO Check if consbonds are complementary
        if self.PS.getMyDiscreteParameter(tc.SOLVER).value is None:
            self.PS.getMyDiscreteParameter(tc.SOLVER).value = random.randint(0, self.N - 1)
        conf_idx = self.PS.getMyDiscreteParameter(tc.SOLVER).value
        self.PS.getMyDiscreteParameter(tc.SOLVER).solutionCounter.recordState()

        if self.config["IK_FlappingSolver"].get("GenerationMode") == "intime":
            for atom_idx in range(len(self.nodes)):
                self.conformer[atom_idx][:] = self.starting_geom[atom_idx][:]
            for idx_frag, c_frag in enumerate(self.appropriate_pairs):
                if conf_idx % 2 ** (idx_frag + 1) >= 2 ** idx_frag:
                    at0 = self.conformer[self.nodes.index(c_frag[0])]
                    at1 = self.conformer[self.nodes.index(c_frag[1])]
                    at2 = self.conformer[self.nodes.index(c_frag[len(c_frag) - 2])]
                    at3 = self.conformer[self.nodes.index(c_frag[len(c_frag) - 1])]

                    av = at2 - at1
                    av /= norm(av)
                    bv = at0 - at1
                    cv = np.cross(av, bv)
                    cv /= norm(cv)
                    bv = np.cross(cv, av)

                    c_frame = np.array([av, bv, cv])
                    tordiff = IK_Math.gettorsion_vec(bv, av, at3 - at2)
                    corr_matrix = np.array([[1, 0, 0],
                                            [0, cos(tordiff / 2), sin(tordiff / 2)],
                                            [0, -sin(tordiff / 2), cos(tordiff / 2)]])
                    corr_matrix = inv(c_frame) @ corr_matrix @ c_frame
                    bv = corr_matrix @ bv
                    cv = corr_matrix @ cv
                    c_frame = np.array([av, bv, cv])
                    reflection_matrix = np.array([[1, 0, 0],
                                                  [0, 1, 0],
                                                  [0, 0, -1]])
                    reflection_matrix = inv(c_frame) @ reflection_matrix @ c_frame

                    for atom in c_frag[2:len(c_frag) - 2]:
                        self.conformer[self.nodes.index(atom)][:] = at1 + reflection_matrix @ \
                                                                       (self.conformer[self.nodes.index(atom)] - at1)

        for item in self.PS:
            if item.isDependent() and item.isContinuous():
                idx = [item.sides[0], item.atoms[0], item.atoms[1], item.sides[1]]
                for i in range(len(idx)):
                    idx[i] = self.nodes.index(idx[i])
                if self.config["IK_FlappingSolver"].get("GenerationMode") == "intime":
                    logger.debug("Required = %f; Actual = %f" % (item.value,
                                                             IK_Math.gettorsion([self.conformer[idx[0]],
                                                                                 self.conformer[idx[1]],
                                                                                 self.conformer[idx[2]],
                                                                                 self.conformer[idx[3]]])))
                else:
                    logger.debug("Required = %f; Actual = %f" % (item.value,
                                                                 IK_Math.gettorsion([self.conformers[conf_idx][idx[0]],
                                                                                     self.conformers[conf_idx][idx[1]],
                                                                                     self.conformers[conf_idx][idx[2]],
                                                                                     self.conformers[conf_idx][idx[3]]])))
        return IK_ParameterSet()
