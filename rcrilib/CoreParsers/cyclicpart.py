from copy import copy, deepcopy
import networkx as nx
import numpy as np
from numpy.linalg import norm, inv

from rcrilib.Helpers import IK_ParameterSet, getfloat, getbool, cause, IK_GeomValidator, geomcheckTree_cycpart, \
    createLogger, IK_Math, IK_CompoundCounter
from .problem import IK_Problem
from .linkingbond import IK_LinkingBond
from .geomunit import IK_GeomUnit
from rcrilib.Helpers import target_class as tc

logger = createLogger("IK_CyclicPart")

class IK_CyclicPart:
    def __init__(self, mygraph, fullmolgraph, atoms_xyz, bonds_num, bonds_type, config):
        self.G = deepcopy(mygraph)
        self.config = config
        for node in list(self.G.nodes):
            self.G.nodes[node]['xyz'] = np.array([atoms_xyz[node][0], atoms_xyz[node][1], atoms_xyz[node][2]])

        self.getcyclegraph()  # CREATE self.CG
        logger.debug("Minimum cycle basis (text line <-> set of nodes of G graph <-> node of CG graph):")
        for i in range(len(self.mcb)):
            logger.debug("%d) %s" % (i + 1, repr(self.mcb[i])))

        self.reducecycles(bonds_num, bonds_type)  # CREATE self.IG
        logger.debug("Linking atoms of independent parts (attributes of IG graph's edges):")
        for edge in list(self.IG.edges()):
            logger.debug("  Link %d - %d : %s" % (edge[0], edge[1], repr(self.IG[edge[0]][edge[1]]['link'])))

        self.redirect_ig()  # CREATE self.DG AND IK_Problem OBJECTS
        self.buildPS()
        if getbool("DoValidation", "IK_CyclicPart", self.config):
            self.validator = IK_GeomValidator(self.G, self.FullPS)
        self.geomunit = IK_GeomUnit(self.G, self.DG, self.headnode, self.config, 0, 0)
        self.ccounter = IK_CompoundCounter()
        self.ccounter += self.geomunit.ccounter

    def getcyclegraph(self):
        mcb = [set(c) for c in nx.minimum_cycle_basis(self.G)]
        edges = []
        curcyc = 0
        i = 1
        while curcyc < len(mcb):
            if i == len(mcb):
                curcyc += 1
                i = curcyc + 1
                continue

            if len(mcb[curcyc].intersection(mcb[i])) > 2:
                mcb[curcyc].update(mcb[i])
                del mcb[i]
                curcyc = 0
                i = 1
            else:
                i += 1
        for i in range(len(mcb)):
            for j in range(i + 1, len(mcb)):
                if len(mcb[i].intersection(mcb[j])) > 0:
                    edges.append([i, j])
        self.CG = nx.Graph()
        if len(edges) == 0:
            if len(mcb) == 1:
                self.CG.add_node(0)
            else:
                raise Exception("No edges and not 1 node!?")
        else:
            self.CG.add_edges_from(edges)
        self.mcb = deepcopy(mcb)

    def isrotatable(self, at1, at2, bonds_num, bonds_type):
        for i in range(len(bonds_num)):
            if at1 in bonds_num[i] and at2 in bonds_num[i]:
                if bonds_type[i] == 1:
                    return True
                else:
                    return False

    def reducecycles(self, bonds_num, bonds_type):
        CGcycles = deepcopy(self.CG)
        for br in nx.bridges(self.CG):
            CGcycles.remove_edge(br[0], br[1])
        cg_comp = [set(CGcycles.subgraph(c).nodes()) for c in
                   nx.connected_components(CGcycles)]  # INDEPENDENT NODES OF self.CG
        indepgraphs = []
        indepnodes = []  # NODES OF INDEPENDENT CYCLES OF self.G
        indepedges = []  # EDGES OF INDEPENDENT CYCLES OF self.G
        for comp in cg_comp:
            newitem = set()
            for i in comp:
                newitem.update(self.mcb[i])
            edges = set()
            for edge in list(self.G.edges):
                if edge[0] in newitem and edge[1] in newitem:
                    edges.add(edge)
            newgraph = nx.Graph()
            newgraph.add_edges_from(edges)
            indepnodes.append(set(newgraph.nodes()))
            indepgraphs.append(newgraph)
            indepedges.append(edges)

        igedges = []  # EDGES OF self.IG
        iglinks = []  # ATTRIBUTES OF EDGES OF self.IG
        for i in range(len(cg_comp)):
            for j in range(i + 1, len(cg_comp)):
                if len(indepnodes[i].intersection(indepnodes[j])) > 0:
                    igedges.append([i, j])
                    iglinks.append([[i, j], list(indepnodes[i].intersection(indepnodes[j]))])

        self.IG = nx.Graph()
        if len(igedges) == 0:
            if len(indepgraphs) == 1:
                self.IG.add_node(0)
            else:
                raise Exception("[CyclicPart] No edges but len(indepgraphs) != 1")
        else:
            self.IG.add_edges_from(igedges)

        for item in iglinks:
            self.IG[item[0][0]][item[0][1]]['link'] = item[1]

        for graph in indepgraphs:
            for node in list(graph.nodes()):
                graph.nodes[node]['xyz'] = deepcopy(self.G.nodes[node]['xyz'])
            for edge in list(graph.edges()):
                graph[edge[0]][edge[1]]['type'] = self.isrotatable(edge[0], edge[1], bonds_num, bonds_type)
                graph[edge[0]][edge[1]]['single'] = self.isrotatable(edge[0], edge[1], bonds_num, bonds_type)

        for i in range(len(indepgraphs)):
            self.IG.nodes[i]['graph'] = indepgraphs[i]
            self.IG.nodes[i]['NDOF'] = self.ndof(indepgraphs[i])

    def redirect_ig(self):
        bestDOF = -100
        self.headnode = -1
        for node in list(self.IG.nodes()):
            if bestDOF < self.IG.nodes[node]['NDOF']:
                self.headnode = node
                bestDOF = self.IG.nodes[node]['NDOF']
        logger.debug("Headnode = " + repr(self.headnode))
        logger.debug("Headnode's atoms = " + repr(list(self.IG.nodes[self.headnode]['graph'].nodes())))
        endnodes = []
        for node in list(self.IG):
            if len(list(self.IG.neighbors(node))) == 1:
                endnodes.append(node)
        if self.headnode in endnodes:
            endnodes.remove(self.headnode)

        self.DG = nx.DiGraph(directed=True)
        logger.debug("Endnodes = " + repr(endnodes))
        for enode in endnodes:
            solvepath = nx.shortest_path(self.IG, source=enode, target=self.headnode)
            for i in range(len(solvepath) - 1):
                self.DG.add_edge(solvepath[i + 1], solvepath[i])
        if len(self.IG.nodes) == 1:
            self.DG.add_node(0)

        for node in list(self.DG.nodes()):
            self.DG.nodes[node]['dep_max'] = len(list(self.DG.neighbors(node)))

        for node in list(self.DG.nodes()):
            links = []
            for neigh_node in list(self.DG.neighbors(node)):
                links.append(self.create_linkingbond(node, neigh_node))
            self.DG.nodes[node]['linking_bonds'] = links

        self.DG = self.DG.reverse(copy=False)

        #GOING FROM ENDS TO HEADNODE
        for node in list(self.DG.nodes()):
            self.DG.nodes[node]['dep_cur'] = 0
            self.DG.nodes[node]['ikprob'] = None
        done = False
        onemorerun = True
        while not done:
            for node in list(self.DG.nodes()):
                if self.DG.nodes[node]['dep_cur'] == self.DG.nodes[node]['dep_max']:
                    if self.DG.nodes[node]['ikprob'] is None:
                        self.DG.nodes[node]['ikprob'] = IK_Problem(self.IG.nodes[node]['graph'], self.config,
                                                                   consbonds=self.DG.nodes[node]['linking_bonds'])
                        self.sync_types(node)
                    elif not self.DG.nodes[node]['ikprob'].recheck_method():
                        self.sync_types(node)
                        onemorerun = True

                    if len(list(self.DG.neighbors(node))) == 0:
                        done = True
                        if onemorerun:
                            onemorerun = False
                            done = False
                            for node in list(self.DG.nodes()):
                                self.DG.nodes[node]['dep_cur'] = 0
                        continue
                    nextnode = list(self.DG.neighbors(node))[0]
                    self.DG.nodes[nextnode]['dep_cur'] += 1
                    self.DG.nodes[node]['dep_cur'] += 1 # TO AVOID REPEATED ACTIONS FOR THE SAME NODE


        for node in list(self.DG.nodes()):
            self.DG.nodes[node]['ikprob'].construct_solver()
            self.DG.nodes[node]['graph'] = self.IG.nodes[node]['graph']

        # SET ATTRIBUTES FOR LINKING TORSIONS
        for node in list(self.DG.nodes()):
            for lb in self.DG.nodes[node]['linking_bonds']:
                logger.info("node = %d; indep = %d; dep = %d" % (node, lb.ind_idx, lb.dep_idx))
                ps = self.DG.nodes[lb.ind_idx]['ikprob'].getPS()

    def startDiscreteRun(self):
        self.geomunit.startDiscreteRun()

    def setname(self, name):
        self.myname = name

    def done_all_discr(self):
        return self.geomunit.done_all_discr()

    def resetDDOF(self):
        self.geomunit.resetDDOF()

    def sync_types(self, node):
        for nb in self.IG.neighbors(node):
            link = self.IG[node][nb]['link']
            if len(link) == 2:
                if not self.IG.nodes[node]['graph'][link[0]][link[1]]['type']:
                    self.IG.nodes[nb]['graph'][link[0]][link[1]]['type'] = False

    def create_linkingbond(self, node, nbnode):
        # node=side1 - dependent
        # nbnode=side2 - independent
        depbond = self.IG[node][nbnode]['link']
        nodegraph = self.IG.nodes[node]['graph']
        nbgraph = self.IG.nodes[nbnode]['graph']

        if len(depbond) == 2:
            node_sideat = [list(nodegraph.neighbors(depbond[0])), list(nodegraph.neighbors(depbond[1]))]
            nb_sideat = [list(nbgraph.neighbors(depbond[0])), list(nbgraph.neighbors(depbond[1]))]
            for item in [node_sideat[0], nb_sideat[0]]:
                item.remove(depbond[1])
            for item in [node_sideat[1], nb_sideat[1]]:
                item.remove(depbond[0])
            linkingbond = IK_LinkingBond(depbond, [node_sideat[0][0], node_sideat[1][0]],
                                         [nb_sideat[0][0], nb_sideat[1][0]],
                                         nbnode, node)
        elif len(depbond) == 1:
            node_sideat = list(nodegraph.neighbors(depbond[0]))
            nb_sideat = list(nbgraph.neighbors(depbond[0]))
            linkingbond = IK_LinkingBond(depbond, [node_sideat[0], node_sideat[1]],
                                         [nb_sideat[0], nb_sideat[1]],
                                         nbnode, node)
        else:
            raise Exception("Trying to link >2 atoms")
        linkingbond.readDependence(self.G)
        return linkingbond

    def buildPS(self):
        self.PS = IK_ParameterSet()
        for node in list(self.DG.nodes()):
            self.DG.nodes[node]['PS'] = self.DG.nodes[node]['ikprob'].getPS()
            self.PS += self.DG.nodes[node]['PS']

        self.FullPS = IK_ParameterSet()
        for node in list(self.DG.nodes()):
            self.FullPS += self.DG.nodes[node]['ikprob'].getFullPS()

    def getPS(self):
        return self.PS.getPS(excludeDep=True)

    def getFullPS(self):
        return self.FullPS.getPS()

    def ndof(self, G):
        """
        TODO ACCOUNT FOR NUMBER OF NEIGHBORS IN self.IG
        """
        addconstr = 0
        ncycles = nx.number_of_edges(self.G) - nx.number_of_nodes(self.G) + 1
        nnodes = G.number_of_nodes()
        return nnodes - 3 - 3 * ncycles - addconstr

    def applyPS(self, increase_discrete = False):
        redoPS = self.geomunit.applyPS(increase_discrete = increase_discrete)
        self.ccounter.checkPS(redoPS)

        if not redoPS.success:
            self.PS.cause = redoPS.cause
            return redoPS.getPS(errorcause=redoPS.cause)

        if getbool("DoValidation", "IK_CyclicPart", self.config):
            good = self.validator.validate(self.G)
            logger.info("Validator returned " + repr(good))
            if not good:
                raise Exception("Failed to obtain IK solution")

        logger.info("Linking is finished!!")
        self.PS.cause = cause.success
        return IK_ParameterSet()

    def setconfig(self, fullmolgraph, dfg, myfrag):
        for i, node in enumerate(self.G.nodes()):
            innerneighb = self.geomunit.getneighbors(node)
            self.G.nodes[node]['innerframe'] = innerneighb

            at1 = self.G.nodes[node]['xyz']
            at0 = self.G.nodes[innerneighb[0]]['xyz']
            at2 = self.G.nodes[innerneighb[1]]['xyz']

            xv = at0 - at1
            yv = at2 - at1
            zv = IK_Math.gs_rand(xv, yv)

            self.G.nodes[node]['outerNB'] = []
            self.G.nodes[node]['localframeXYZ'] = []
            for molnode in fullmolgraph.neighbors(node):
                if molnode not in innerneighb: #and not self.G.has_node(molnode):
                    self.G.nodes[node]['outerNB'].append(molnode)
                    bondvec = fullmolgraph.nodes[molnode]['xyz'] - at1
                    self.G.nodes[node]['localframeXYZ'].append(np.array([bondvec @ xv,
                                                                         bondvec @ yv,
                                                                         bondvec @ zv]))
            self.G.nodes[node]['outerXYZ'] = [0] * len(self.G.nodes[node]['outerNB'])

        if len(list(dfg.neighbors(myfrag))) == 1:
            self.recalc_frame = True
            otherfrag = list(dfg.neighbors(myfrag))[0]
            bond = dfg[myfrag][otherfrag]['param'][0].atoms
            sides = dfg[myfrag][otherfrag]['param'][0].sides
            if bond[1] in dfg.nodes[myfrag]['atoms']:
                self.at0_idx = bond[0]
                self.at1_idx = bond[1]
                self.at2_idx = sides[1]
            elif bond[0] in dfg.nodes[myfrag]['atoms']:
                self.at0_idx = bond[1]
                self.at1_idx = bond[0]
                self.at2_idx = sides[0]
        elif len(list(dfg.neighbors(myfrag))) == 0:
            self.recalc_frame = False
        else:
            raise Exception("dfg IS NOT A TREE!!!")

        #EXTEND G WITH FIRST NON-CYCLIC NEIGHBORS
        self.extG = nx.Graph()
        self.extG.add_nodes_from(self.G.nodes)
        self.extG.add_edges_from(self.G.edges)
        for cnode_idx, cnode in enumerate(self.G.nodes()):
            for node in self.G.nodes[cnode]['outerNB']:
                if not self.extG.has_edge(cnode, node):
                    self.extG.add_edge(cnode, node)

        self.geomunit.setExtG(self.extG)

    def get_outer_xyz(self, cnode, Gidx):
        for i, index in enumerate(self.G.nodes[cnode]['outerNB']):
            if index == Gidx:
                return self.G.nodes[cnode]['outerXYZ'][i]
        raise Exception("Something's wrong. I can feel it")

    def updateattributes(self, dfg, myfrag, loc_coord):
        if self.recalc_frame:
            at0 = self.get_outer_xyz(self.at1_idx, self.at0_idx)
            at1 = self.G.nodes[self.at1_idx]['xyz']
            at2 = self.G.nodes[self.at2_idx]['xyz']
            xv = at1 - at0
            yv = at2 - at1
            zv = IK_Math.gs_rand(xv, yv)
        else:
            at1 = np.zeros(3)
            xv = np.array([1, 0, 0])
            yv = np.array([0, 1, 0])
            zv = np.array([0, 0, 1])

        for cnode_idx, cnode in enumerate(self.G.nodes()):
            loc_coord[cnode] = np.array([(self.G.nodes[cnode]['xyz'] - at1) @ xv,
                                         (self.G.nodes[cnode]['xyz'] - at1) @ yv,
                                         (self.G.nodes[cnode]['xyz'] - at1) @ zv,
                                         1])
            for onode_idx, onode in enumerate(self.G.nodes[cnode]['outerNB']):
                if onode in dfg.nodes[myfrag]['atoms'] and onode not in self.G.nodes():
                    loc_coord[onode] = np.array([(self.G.nodes[cnode]['outerXYZ'][onode_idx] - at1) @ xv,
                                                 (self.G.nodes[cnode]['outerXYZ'][onode_idx] - at1) @ yv,
                                                 (self.G.nodes[cnode]['outerXYZ'][onode_idx] - at1) @ zv,
                                                 1])

        startframe = np.zeros((4, 4))
        startframe[:3, 0] = xv
        startframe[:3, 1] = yv
        startframe[:3, 2] = zv
        startframe[:3, 3] = at1
        startframe[3, 3] = 1
        dfg.nodes[myfrag]['frame'] = startframe
        for nb in list(dfg.neighbors(myfrag)):
            endbond = dfg[myfrag][nb]['param'][0].atoms
            endsides = dfg[myfrag][nb]['param'][0].sides

            if endbond[1] in dfg.nodes[nb]['atoms']:
                at0 = self.get_outer_xyz(endbond[0], endbond[1])
                at1 = self.G.nodes[endbond[0]]['xyz']
                at2 = self.G.nodes[endsides[0]]['xyz']
            elif endbond[0] in dfg.nodes[nb]['atoms']:
                at0 = self.get_outer_xyz(endbond[1], endbond[0])
                at1 = self.G.nodes[endbond[1]]['xyz']
                at2 = self.G.nodes[endsides[1]]['xyz']

            xv = at0 - at1
            yv = at2 - at1
            zv = IK_Math.gs_rand(xv, yv)

            newframe = np.zeros((4, 4))
            newframe[:3, 0] = xv
            newframe[:3, 1] = yv
            newframe[:3, 2] = zv
            newframe[:3, 3] = at0
            newframe[3, 3] = 1
            dfg[myfrag][nb]['transition'] = inv(startframe) @ newframe

    def perturb_geometry(self):
        self.geomunit.perturb_geometry()

    def check_overlap(self):
        return geomcheckTree_cycpart(self.G, self.extG, getfloat("MinDistance", "IK_CyclicPart", self.config))
