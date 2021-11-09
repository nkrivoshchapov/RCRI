from copy import copy, deepcopy
import networkx as nx
import numpy as np
from numpy.linalg import norm, inv

from rcrilib.Helpers import IK_ParameterSet, getfloat, getbool, cause, IK_GeomValidator, geomcheckTree_cycpart, \
    createLogger, IK_Math, IK_CompoundCounter
from .problem import IK_Problem, ikalg
from .linkingbond import IK_LinkingBond
from .geomunit import IK_GeomUnit
from .fakebond import IK_FakeBond

logger = createLogger("IK_CyclicPart")

class IK_CyclicPart:
    def __init__(self, mygraph, fullmolgraph, atoms_xyz, bonds_num, bonds_type, config):
        self.G = deepcopy(mygraph)
        self.config = config
        for node in list(self.G.nodes):
            self.G.nodes[node]['xyz'] = np.array([atoms_xyz[node][0], atoms_xyz[node][1], atoms_xyz[node][2]])

        self.getcyclegraph(bonds_num, bonds_type)  # CREATE self.IG
        logger.debug("Minimum cycle basis (text line <-> set of nodes of G graph <-> node of CG graph):")
        for i in range(len(self.mcb)):
            logger.debug("%d) %s" % (i + 1, repr(self.mcb[i])))

        # self.reducecycles(bonds_num, bonds_type)
        logger.debug("Linking atoms of independent parts (attributes of IG graph's edges):")
        for edge in list(self.IG.edges()):
            logger.debug("  Link %d - %d : %s" % (edge[0], edge[1], repr(self.IG[edge[0]][edge[1]]['link'])))

        self.redirect_ig() # CREATE self.DG AND IK_Problem OBJECTS
        self.buildPS()
        if getbool("DoValidation", "IK_CyclicPart", self.config) or getbool("DoValidation", "IK_GeomUnit", self.config):
            if getbool("AllowAnglePerturbation", "IK_TLCSolver", self.config):
                for graph_idx in self.IG.nodes():
                    curgraph = self.IG.nodes[graph_idx]['graph']
                    for node in curgraph.nodes():
                        if "shared_vangle" in curgraph.nodes[node]:
                            self.G.nodes[node]['shared_vangle'] = curgraph.nodes[node]['shared_vangle']
            if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config):
                for graph_idx in self.IG.nodes():
                    curgraph = self.IG.nodes[graph_idx]['graph']
                    for edge in curgraph.edges():
                        if "shared_length" in curgraph[edge[0]][edge[1]]:
                            self.G[edge[0]][edge[1]]['shared_length'] = curgraph[edge[0]][edge[1]]['shared_length']

        if getbool("DoValidation", "IK_CyclicPart", self.config):
            self.validator = IK_GeomValidator(self.G, self.FullPS)
        self.geomunit = IK_GeomUnit(self.G, self.DG, self.headnode, self.config, 0, 0)
        self.ccounter = IK_CompoundCounter()
        self.ccounter += self.geomunit.ccounter

    def getcyclegraph(self, bonds_num, bonds_type):
        first_mcb = [set(c) for c in nx.minimum_cycle_basis(self.G)]
        tempG = nx.Graph()
        for i in range(len(first_mcb)):
            tempG.add_node(i)
            tempG.nodes[i]['atomset'] = first_mcb[i]
            if len(first_mcb[i]) <= 5:
                tempG.nodes[i]['ndof'] = -1 # THEY MUST BE UNITED
            else:
                tempG.nodes[i]['ndof'] = len(first_mcb[i]) - 6

        for i in range(tempG.number_of_nodes()):
            for j in range(i + 1, tempG.number_of_nodes()):
                if len(tempG.nodes[i]['atomset'].intersection(tempG.nodes[j]['atomset'])) > 0:
                    tempG.add_edge(i, j)
        clean_run = False
        while not clean_run:
            clean_run = True
            for node in list(tempG.nodes()):
                if tempG.nodes[node]['ndof'] == -1:
                    for nbnode in list(tempG.neighbors(node)):
                        spiro_link = (len(tempG.nodes[node]['atomset'].
                                          intersection(tempG.nodes[nbnode]['atomset'])) == 1)
                        lack_dofs = (len(tempG.nodes[nbnode]['atomset']) - 6 - (len(tempG.nodes[node]['atomset'].
                                          intersection(tempG.nodes[nbnode]['atomset'])) - 1) < 0)
                        logger.info("node atoms = %s\nnbnode atoms = %s\nspiro criterion = %d(%s)\n"
                                     "lack criterion = %d(%s)" % (
                            repr(tempG.nodes[node]['atomset']),
                            repr(tempG.nodes[nbnode]['atomset']),
                            len(tempG.nodes[node]['atomset'].intersection(tempG.nodes[nbnode]['atomset'])),
                            repr(spiro_link),
                            len(tempG.nodes[nbnode]['atomset']) - 6 - (len(tempG.nodes[node]['atomset'].
                                                   intersection(tempG.nodes[nbnode]['atomset'])) - 1), repr(lack_dofs)
                        ))
                        if not spiro_link and lack_dofs:
                            # Should unite 'node' with 'nbnode'
                            tempG.nodes[node]['atomset'].update(tempG.nodes[nbnode]['atomset'])
                            for other_nb in tempG.neighbors(nbnode):
                                if other_nb != node:
                                    tempG.add_edge(other_nb, node)
                            tempG.remove_node(nbnode)
                            clean_run = False
                    if not clean_run:
                        break

        self.mcb = []
        for node in tempG.nodes():
            self.mcb.append(tempG.nodes[node]['atomset'])
        self.IG = nx.Graph()
        for i in range(len(self.mcb)):
            for j in range(i + 1, len(self.mcb)):
                common_atoms = self.mcb[i].intersection(self.mcb[j])
                if len(common_atoms) > 0:
                    self.IG.add_edge(i, j, link=list(common_atoms))
        if self.IG.number_of_edges() == 0 and len(self.mcb) == 1:
            self.IG.add_node(0)
        assert self.IG.number_of_nodes() != 0, "Cycle graph has 0 nodes"

        for i in range(len(self.mcb)):
            fraggraph = nx.Graph()
            for edge in list(self.G.edges):
                if edge[0] in self.mcb[i] and edge[1] in self.mcb[i]:
                    issingle = self.isrotatable(edge[0], edge[1], bonds_num, bonds_type)
                    fraggraph.add_edge(*edge, type=issingle, single=issingle)
            for node in fraggraph.nodes():
                fraggraph.nodes[node]['xyz'] = deepcopy(self.G.nodes[node]['xyz'])
            self.IG.nodes[i]['graph'] = fraggraph
            self.IG.nodes[i]['NDOF'] = self.ndof(fraggraph)
        # TODO Detect topology of a sphere (3D cage)

        IGcycles = nx.Graph()
        IGcycles.add_edges_from(self.IG.edges)
        for br in nx.bridges(self.IG):
            IGcycles.remove_edge(br[0], br[1])
        self.preDG = nx.Graph()
        if len(self.IG.edges) == 0:
            self.preDG.add_nodes_from(self.IG.nodes)
        else:
            self.preDG.add_edges_from(self.IG.edges)
        ringgroups = [self.preDG.subgraph(IGcycles.subgraph(c).nodes()) for c in nx.connected_components(IGcycles)]

        for cursubgraph in ringgroups:
            if nx.number_of_nodes(cursubgraph) > 1: # Then it is a ring in of condensed cycles
                cyclicedges = set(cursubgraph.edges)
                while nx.number_of_nodes(cursubgraph) - nx.number_of_edges(cursubgraph) != 1:
                    cyclicedges.difference_update(set(nx.bridges(cursubgraph)))
                    # logger.info("Bridges = " + repr(set(nx.bridges(cursubgraph))))
                    # logger.info("cyclicedges = " + repr(cyclicedges))
                    mindof = None
                    minnode = None
                    for node in cursubgraph.nodes():
                        ndof = self.IG.nodes[node]['NDOF']
                        if mindof is None or mindof > ndof:
                            for edge in cyclicedges:
                                if node in edge:
                                    mindof = ndof
                                    minnode = node
                                    break
                    nbs = list(cursubgraph.neighbors(minnode))
                    maxdof = None
                    maxnode = None
                    for nb in nbs:
                        ndof = self.IG.nodes[nb]['NDOF']
                        if maxdof is None or maxdof < ndof:
                            maxdof = ndof
                            maxnode = nb
                    self.preDG.remove_edge(maxnode, minnode)
                    # logger.error("Removed edge %d - %d" % (maxnode, minnode))
                    if (maxnode, minnode) in cyclicedges:
                        cyclicedges.remove((maxnode, minnode))
                    else:
                        cyclicedges.remove((minnode, maxnode))

    def isrotatable(self, at1, at2, bonds_num, bonds_type):
        for i in range(len(bonds_num)):
            if at1 in bonds_num[i] and at2 in bonds_num[i]:
                if bonds_type[i] == 1:
                    return True
                else:
                    return False

    # def reducecycles(self, bonds_num, bonds_type):
    #     CGcycles = deepcopy(self.CG)
    #     for br in nx.bridges(self.CG):
    #         CGcycles.remove_edge(br[0], br[1])
    #     cg_comp = [set(CGcycles.subgraph(c).nodes()) for c in
    #                nx.connected_components(CGcycles)]  # INDEPENDENT NODES OF self.CG
    #     indepgraphs = []
    #     indepnodes = []  # NODES OF INDEPENDENT CYCLES OF self.G
    #     indepedges = []  # EDGES OF INDEPENDENT CYCLES OF self.G
    #     for comp in cg_comp:
    #         newitem = set()
    #         for i in comp:
    #             newitem.update(self.mcb[i])
    #         edges = set()
    #         for edge in list(self.G.edges):
    #             if edge[0] in newitem and edge[1] in newitem:
    #                 edges.add(edge)
    #         newgraph = nx.Graph()
    #         newgraph.add_edges_from(edges)
    #         indepnodes.append(set(newgraph.nodes()))
    #         indepgraphs.append(newgraph)
    #         indepedges.append(edges)
    #
    #     igedges = []  # EDGES OF self.IG
    #     iglinks = []  # ATTRIBUTES OF EDGES OF self.IG
    #     for i in range(len(cg_comp)):
    #         for j in range(i + 1, len(cg_comp)):
    #             if len(indepnodes[i].intersection(indepnodes[j])) > 0:
    #                 igedges.append([i, j])
    #                 iglinks.append([[i, j], list(indepnodes[i].intersection(indepnodes[j]))])
    #
    #     self.IG = nx.Graph()
    #     if len(igedges) == 0:
    #         if len(indepgraphs) == 1:
    #             self.IG.add_node(0)
    #         else:
    #             raise Exception("[CyclicPart] No edges but len(indepgraphs) != 1")
    #     else:
    #         self.IG.add_edges_from(igedges)
    #
    #     for item in iglinks:
    #         self.IG[item[0][0]][item[0][1]]['link'] = item[1]
    #
    #     for graph in indepgraphs:
    #         for node in list(graph.nodes()):
    #             graph.nodes[node]['xyz'] = deepcopy(self.G.nodes[node]['xyz'])
    #         for edge in list(graph.edges()):
    #             graph[edge[0]][edge[1]]['type'] = self.isrotatable(edge[0], edge[1], bonds_num, bonds_type)
    #             graph[edge[0]][edge[1]]['single'] = self.isrotatable(edge[0], edge[1], bonds_num, bonds_type)
    #
    #     for i in range(len(indepgraphs)):
    #         self.IG.nodes[i]['graph'] = indepgraphs[i]
    #         self.IG.nodes[i]['NDOF'] = self.ndof(indepgraphs[i])

    def redirect_ig(self):
        bestDOF = None
        self.headnode = None
        for node in list(self.IG.nodes()):
            if bestDOF is None or bestDOF < self.IG.nodes[node]['NDOF']:
                self.headnode = node
                bestDOF = self.IG.nodes[node]['NDOF']
        logger.info("Headnode = " + repr(self.headnode))
        logger.info("Headnode's atoms = " + repr(list(self.IG.nodes[self.headnode]['graph'].nodes())))
        endnodes = []
        for node in list(self.preDG):
            if len(list(self.preDG.neighbors(node))) == 1:
                endnodes.append(node)
        if self.headnode in endnodes:
            endnodes.remove(self.headnode)

        self.DG = nx.DiGraph(directed=True)
        logger.debug("Endnodes = " + repr(endnodes))

        for enode in endnodes:
            solvepath = nx.shortest_path(self.preDG, source=enode, target=self.headnode)
            for i in range(len(solvepath) - 1):
                if not self.DG.has_edge(solvepath[i + 1], solvepath[i]):
                    self.DG.add_edge(solvepath[i + 1], solvepath[i])
                    self.DG.nodes[solvepath[i]]['dg_order'] = len(solvepath) - i - 1
                    logger.info("Adding edge to DG : %d -> %d" % (solvepath[i + 1], solvepath[i]))
                    logger.info("Node %d has order %d" % (solvepath[i], self.DG.nodes[solvepath[i]]['dg_order']))
        if len(self.preDG.nodes) == 1:
            self.DG.add_node(0)
        self.DG.nodes[self.headnode]['dg_order'] = 0

        for node in self.DG.nodes:
            self.DG.nodes[node]['dep_max'] = len(list(self.DG.neighbors(node)))

        self.LG = nx.DiGraph(directed=True)
        self.LG.add_edges_from(self.DG.edges)
        if len(self.DG.edges) == 0:
            self.LG.add_nodes_from(self.DG.nodes)

        for edge in self.IG.edges:
            if not self.LG.has_edge(edge[0], edge[1]) and not self.LG.has_edge(edge[1], edge[0]):
                if nx.has_path(self.DG, source=edge[0], target=edge[1]):
                    self.LG.add_edge(edge[0], edge[1])
                elif nx.has_path(self.DG, source=edge[1], target=edge[0]):
                    self.LG.add_edge(edge[1], edge[0])
                else:
                    newedge = None
                    if self.DG.nodes[edge[0]]['dg_order'] > self.DG.nodes[edge[1]]['dg_order']:
                        newedge = [edge[1], edge[0]]
                    elif self.DG.nodes[edge[0]]['dg_order'] < self.DG.nodes[edge[1]]['dg_order']:
                        newedge = [edge[0], edge[1]]
                    else:
                        newedge = [edge[0], edge[1]] # These choices will not really matter
                    self.LG.add_edge(newedge[0], newedge[1])
                    logger.info("Hacking the LG path %d - %d" % (newedge[0], newedge[1]))
        logger.info("DG Edges = " + repr(self.DG.edges))
        logger.info("LG Edges = " + repr(self.LG.edges))

        for node in self.LG.nodes:
            links = []
            for neigh_node in list(self.LG.neighbors(node)):
                links.append(self.create_linkingbond(node, neigh_node))
            self.LG.nodes[node]['linking_bonds'] = links

        self.DG = self.DG.reverse(copy=False)
        self.LG = self.LG.reverse(copy=False)

        #GOING FROM ENDS TO HEADNODE
        for node in self.DG.nodes:
            self.DG.nodes[node]['dep_cur'] = 0
            self.DG.nodes[node]['ikprob'] = None
        done = False
        onemorerun = True
        while not done:
            for node in list(self.DG.nodes()):
                if self.DG.nodes[node]['dep_cur'] == self.DG.nodes[node]['dep_max']:
                    if self.DG.nodes[node]['ikprob'] is None:
                        self.DG.nodes[node]['ikprob'] = IK_Problem(self.IG.nodes[node]['graph'], self.config,
                                                                   consbonds=self.LG.nodes[node]['linking_bonds'])
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
        # TODO If we got two neighboring IK_FlappingSolvers (or IK_45Solvers) give warning that code is not optimized

        # ASSIGN SHARED VARIABLES TO BOND LENGTHS AND VALENCE ANGLES
        if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config) or \
           getbool("AllowAnglePerturbation", "IK_TLCSolver", self.config):
            if getbool("DoValidation", "IK_CyclicPart", self.config) or \
               getbool("DoValidation", "IK_Molecule", self.config) or \
               getbool("DoValidation", "IK_Problem", self.config):
                for node in self.DG.nodes():
                    if self.DG.nodes[node]['ikprob'].method == ikalg.TLC:
                        if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config):
                            for edge in self.IG.nodes[node]['graph'].edges():
                                svalue = self.DG.nodes[node]['ikprob'].share_length(edge)

            for depnode in self.DG.nodes():
                lbs = self.LG.nodes[depnode]['linking_bonds']
                for lb in lbs:
                    indnode = lb.ind_idx
                    if lb.isFake():
                        if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config):
                            for bond in lb.sbonds + lb.ibonds:
                                bond_gidx = [bond['atoms'][0]['G_idx'], bond['atoms'][1]['G_idx']]
                                svalue = self.DG.nodes[indnode]['ikprob'].share_length(bond_gidx)
                                self.DG.nodes[depnode]['ikprob'].assign_shared_length(bond_gidx, svalue)
                        if getbool("AllowAnglePerturbation", "IK_TLCSolver", self.config):
                            for atom in lb.iatoms:
                                gidx = atom['G_idx']
                                svalue = self.DG.nodes[indnode]['ikprob'].share_vangle(gidx)
                                self.DG.nodes[depnode]['ikprob'].assign_shared_vangle(gidx, svalue)
                    elif getbool("AllowBondPerturbation", "IK_TLCSolver", self.config) and len(lb.bond) == 2:
                        svalue = self.DG.nodes[indnode]['ikprob'].share_length(lb.bond)
                        self.DG.nodes[depnode]['ikprob'].assign_shared_length(lb.bond, svalue)

        # INITIALIZE SOLVER OBJECTS
        for node in list(self.DG.nodes()):
            self.DG.nodes[node]['graph'] = self.IG.nodes[node]['graph']
            self.DG.nodes[node]['ikprob'].construct_solver()

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
            elif len(link) > 2:
                for edge in self.IG.nodes[node]['graph'].edges():
                    if edge[0] in link and edge[1] in link and \
                            not self.IG.nodes[node]['graph'][edge[0]][edge[1]]['type']:
                        self.IG.nodes[nb]['graph'][edge[0]][edge[1]]['type'] = False

    def create_linkingbond(self, node, nbnode):
        # node=side1 - dependent
        # nbnode=side2 - independent
        depbond = self.IG[node][nbnode]['link'] # List of linking atoms
        nodegraph = self.IG.nodes[node]['graph']
        nbgraph = self.IG.nodes[nbnode]['graph']

        if len(depbond) == 2:
            node_sideat = [list(nodegraph.neighbors(depbond[0])), list(nodegraph.neighbors(depbond[1]))]
            nb_sideat = [list(nbgraph.neighbors(depbond[0])), list(nbgraph.neighbors(depbond[1]))]
            for item in [node_sideat[0], nb_sideat[0]]:
                item.remove(depbond[1])
            for item in [node_sideat[1], nb_sideat[1]]:
                item.remove(depbond[0])
            lb = IK_LinkingBond(depbond, [node_sideat[0][0], node_sideat[1][0]],
                                         [nb_sideat[0][0], nb_sideat[1][0]],
                                         nbnode, node)
            lb.readDependence(self.G)
            return lb
        elif len(depbond) == 1:
            node_sideat = list(nodegraph.neighbors(depbond[0]))
            nb_sideat = list(nbgraph.neighbors(depbond[0]))
            lb = IK_LinkingBond(depbond, [node_sideat[0], node_sideat[1]],
                                         [nb_sideat[0], nb_sideat[1]],
                                         nbnode, node)
            lb.readDependence(self.G)
            return lb
        else:
            fb = IK_FakeBond(nodegraph, nbgraph, depbond, nbnode, node, self.config)
            fb.readDependence(self.G)
            return fb

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
        ncycles = nx.number_of_edges(G) - nx.number_of_nodes(G) + 1
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
        self.geomunit.init_polyhedra(fullmolgraph)

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
            for node in self.G.nodes[cnode]['main_outerNB']:
                if not self.extG.has_edge(cnode, node):
                    self.extG.add_edge(cnode, node)
        self.geomunit.setExtG(self.extG)

    def get_outer_xyz(self, cnode, Gidx):
        for i, index in enumerate(self.G.nodes[cnode]['main_outerNB']):
            if index == Gidx:
                return self.G.nodes[cnode]['main_outerXYZ'][i]
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
            for onode_idx, onode in enumerate(self.G.nodes[cnode]['main_outerNB']):
                if onode in dfg.nodes[myfrag]['atoms'] and onode not in self.G.nodes():
                    loc_coord[onode] = np.array([(self.G.nodes[cnode]['main_outerXYZ'][onode_idx] - at1) @ xv,
                                                 (self.G.nodes[cnode]['main_outerXYZ'][onode_idx] - at1) @ yv,
                                                 (self.G.nodes[cnode]['main_outerXYZ'][onode_idx] - at1) @ zv,
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
