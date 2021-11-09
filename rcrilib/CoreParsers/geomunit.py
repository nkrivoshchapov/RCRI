import networkx as nx
from copy import deepcopy
import numpy as np
import builtins

from rcrilib.Helpers import createLogger, IK_ParameterSet, getfloat, getbool, cause, IK_GeomValidator, \
    geomcheckTree_graph, IK_Math, IK_CompoundCounter
from rcrilib.Helpers  import target_class as tc

logger = createLogger("IK_GeomUnit")

class IK_GeomUnit:
    def __init__(self, G, DG, headnode, config, ord, idx, master_unit=None):
        self.G = G
        self.DG = DG
        self.headnode = headnode
        self.config = config
        self.idx = idx
        self.ord = ord
        self.masterunit = master_unit

        tempDG = nx.Graph(directed=True)
        tempDG.add_nodes_from(DG)
        tempDG.add_edges_from(DG.edges)
        sub_heads = set(tempDG.neighbors(headnode))
        tempDG.remove_node(headnode)
        dg_comps = []
        for c in nx.connected_components(tempDG):
            sub_head = (sub_heads & set(c)).pop()
            dg_comps.append([DG.subgraph(c), sub_head])

        self.dep_indicies = []
        self.deps = []
        for i in range(len(dg_comps)):
            dg_sub = dg_comps[i][0]
            sub_head = dg_comps[i][1]
            g_sub_nodes = []
            for node in dg_sub:
                g_sub_nodes += list(DG.nodes[node]['graph'].nodes())
            g_sub = G.subgraph(g_sub_nodes)

            newunit = IK_GeomUnit(g_sub, dg_sub, sub_head, self.config, self.ord+1, i, master_unit=self)
            self.deps.append(newunit)
            self.dep_indicies.append(newunit.headnode)

        self.FullPS = IK_ParameterSet()
        for node in list(self.DG.nodes()):
            self.FullPS += self.DG.nodes[node]['ikprob'].getFullPS()
        if getbool("DoValidation", "IK_GeomUnit", self.config):
            self.validator = IK_GeomValidator(self.G, self.FullPS)

        self.scounter = self.DG.nodes[self.headnode]['ikprob'].scounter
        self.ccounter = IK_CompoundCounter()
        self.ccounter += self.scounter
        for dep in self.deps:
            self.ccounter += dep.ccounter
        self.DG.nodes[self.headnode]['geomunit'] = self
        self.good_ddof = -1

        # logger.error("GU order=%d; idx=%d; nodes = %s" % (self.ord, self.idx, repr(list(self.G.nodes))))
        if self.ord == 0:
            for dep in self.deps:
                dep.finalize_init()
            self.finalize_init()

    def finalize_init(self):
        self.lbs = []
        for lb in self.DG.nodes[self.headnode]['ikprob'].consbonds:
            newlb = {'lb': lb, 'do_linking': (lb.ind_idx in self.dep_indicies)}
            newlb['local_graph'], newlb['global_graph'] = self.transmit_graphs(lb.ind_idx, trymaster=True)
            self.lbs.append(newlb)
        for dep in self.deps:
            dep.finalize_init()

    def init_polyhedra(self, fullgraph):
        if self.ord == 0:
            for node in self.G.nodes:
                self.G.nodes[node]['frames'] = []

        for dep in self.deps:
            dep.init_polyhedra(fullgraph)

        for i, node in enumerate(self.G.nodes):
            createframe = True
            myframeidx = None
            for frame_idx, frame in enumerate(self.G.nodes[node]['frames']):
                accepted = True
                for atom in frame['innerframe']:
                    if not self.G.has_node(atom):
                        accepted = False
                        break
                if accepted:
                    createframe = False
                    myframeidx = frame_idx
                    break

            if createframe:
                innerneighb = list(self.G.neighbors(node))[:2]
                myframeidx = len(self.G.nodes[node]['frames'])
                myframe = {'innerframe': innerneighb, 'outerNB': [], 'localframeXYZ': [], 'users': []}

                at1 = self.G.nodes[node]['xyz']
                at0 = self.G.nodes[myframe['innerframe'][0]]['xyz']
                at2 = self.G.nodes[myframe['innerframe'][1]]['xyz']
                xv = at0 - at1
                yv = at2 - at1
                zv = IK_Math.gs_rand(xv, yv)

                for molnode in fullgraph.neighbors(node):
                    if molnode not in innerneighb:
                        myframe['outerNB'].append(molnode)
                        bondvec = fullgraph.nodes[molnode]['xyz'] - at1
                        myframe['localframeXYZ'].append(np.array([bondvec @ xv,
                                                                  bondvec @ yv,
                                                                  bondvec @ zv]))
                myframe['outerXYZ'] = [None] * len(myframe['outerNB'])
                self.G.nodes[node]['frames'].append(myframe)

            assert myframeidx is not None
            self.G.nodes[node]['frames'][myframeidx]['users'].append(self.headnode)

        if self.ord == 0:
            for cnode_idx, cnode in enumerate(self.G.nodes):
                main_frame_idx = None
                for frameidx, frame in enumerate(self.G.nodes[cnode]['frames']):
                    if self.headnode in frame['users']:
                        main_frame_idx = frameidx
                        break
                assert main_frame_idx is not None
                self.G.nodes[cnode]['main_outerNB'] = self.G.nodes[cnode]['frames'][main_frame_idx]['outerNB']
                self.G.nodes[cnode]['main_outerXYZ'] = self.G.nodes[cnode]['frames'][main_frame_idx]['outerXYZ']

    def __contains__(self, node):
        return self.G.has_node(node)

    def getneighbors(self, node):
        for dep in self.deps:
            if node in dep:
                return dep.getneighbors(node)
        return list(self.G.neighbors(node))[:2]

    def setExtG(self, extG):
        self.extG = extG
        for i in range(len(self.deps)):
            extg_sub_nodes = []
            for node in list(self.deps[i].G.nodes()):
                extg_sub_nodes.append(node)
                for nb in list(extG.neighbors(node)):
                    if nb not in extg_sub_nodes:
                        extg_sub_nodes.append(nb)
            extg_sub = extG.subgraph(extg_sub_nodes)
            self.deps[i].setExtG(extg_sub)

    def startDiscreteRun(self):
        self.ccounter.startDiscreteRun()

    def done_self_discr(self):
        return self.scounter.done_discr()

    def done_all_discr(self):
        if not self.done_self_discr():
            logger.info("order=%d; idx=%d Is stopping from termination atoms=%s" % (self.ord, self.idx, repr(list(self.G.nodes()))))
        return self.ccounter.done_discr()

    def resetDDOF(self):
        logger.info("Excel: line 129 resetDDOF")
        self.ccounter.resetDDOF()

    def applyPS(self, increase_discrete = False):
        # if builtins.loglinecount == 30 and self.ord==1 and self.idx==0:
        #     print("HERE")
        # self.double_check()
        # if builtins.loglinecount == 30 and self.ord == 3 and self.idx==0:
        #     self.masterunit.masterunit.double_check(hard=True)

        if not increase_discrete and getbool("SmartAPS", "IK_GeomUnit", self.config) and not self.FullPS.isChanged():
            # if self.ord==1 and self.idx==0 and builtins.loglinecount == 30:
            #     self.double_check(hard=True)
            logger.info("Leaving GU order=%d; idx=%d unchanged. LS = "%(self.ord, self.idx) + repr(self.FullPS.cause))
            if hasattr(self.scounter, 'discValues'):
                logger.info("Skipping %s order=%d; idx=%d My list = %s cursolution = %s" % (
                                repr(increase_discrete), self.ord, self.idx, repr(self.scounter.discValues),
                                repr(self.scounter.ddof.getValue())))
                # for dep in self.deps:
                #     dep.applyPS()
            # if self.ord==1 and self.idx==0:
            #     self.double_check(hard=True)
            if self.FullPS.cause == cause.success:
                return IK_ParameterSet()
            else:
                ps = self.DG.nodes[self.headnode]['ikprob'].getPS(errorcause=self.FullPS.cause)
                if self.scounter.isDepChanged(ps) or self.FullPS.cause == cause.geomoverlap:
                    self.scounter.resetDDOF()
                return ps
        logger.info("Entering GU order=%d; idx=%d" % (self.ord, self.idx))

        incremented = False
        if increase_discrete:
            if not self.scounter.done_discr():
                self.scounter.nextValue()
                for dep in self.deps:
                    ps = dep.applyPS()
                    if not ps.success:
                        if self.scounter.isDepChanged(ps):
                            self.scounter.resetDDOF()
                        logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
                        self.FullPS.cause = cause.success  # NO WAY TO RETURN A RIGHT PS
                        return ps.getPS(errorcause=ps.cause)
            else:
                curdep = -1
                for i in range(len(self.deps)):
                    if not incremented and not self.deps[i].done_all_discr():
                        curdep = i
                        break
                if curdep == -1:
                    raise Exception("Did not expect -1")
                ps = self.deps[curdep].applyPS(increase_discrete=True)
                logger.info("Excel: line 167 resetDDOF")
                self.scounter.resetDDOF()
                if not ps.success:
                    self.FullPS.cause = cause.success  # NO WAY TO RETURN A RIGHT PS
                    logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
                    return ps.getPS(errorcause=ps.cause)
                incremented = True

                for i in range(len(self.deps)):
                    if i != curdep:
                        ps = self.deps[i].applyPS()
                        if not ps.success:
                            if self.scounter.isDepChanged(ps):
                                self.scounter.resetDDOF()
                            self.FullPS.cause = cause.success  # NO WAY TO RETURN A RIGHT PS
                            logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
                            return ps.getPS(errorcause=ps.cause)
        else:
            for dep in self.deps:
                ps = dep.applyPS()
                if not ps.success:
                    if self.scounter.isDepChanged(ps):
                        self.scounter.resetDDOF()
                    self.FullPS.cause = cause.success # NO WAY TO RETURN A RIGHT PS
                    logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
                    return ps.getPS(errorcause=ps.cause)
        for lb in self.lbs:
            lb['lb'].resolveDependence(lb['local_graph']) # YDACHI BRATAN
        # if self.ord==1 and self.idx==0 and builtins.loglinecount == 27:
        #     self.deps[0].double_check(hard=True)
        #     print("HERE")
        logger.info("Running APS for GU order=%d; idx=%d" % (self.ord, self.idx))
        if hasattr(self.scounter, 'discValues'):
            logger.info("%s order=%d; idx=%d My list = %s cursolution = %s" % (repr(increase_discrete),
                                                                               self.ord, self.idx,
                                                                               repr(self.scounter.discValues),
                                                                               repr(self.scounter.ddof.getValue())))

        redoPS = self.DG.nodes[self.headnode]['ikprob'].applyPS()
        if getbool("DoValidation", "IK_CyclicPart", self.config):
            for dep in self.deps:
                dep.double_check(hard=True)
        # if self.ord == 1 and self.idx == 0 and builtins.loglinecount == 27:

        logger.debug("redoPS = " + repr(redoPS))
        if not redoPS.success:
            self.FullPS.cause = redoPS.cause
            logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
            return redoPS.getPS(errorcause=redoPS.cause) # TODO Don't touch torsions on linking bonds
        if incremented:
            self.scounter.startDiscreteRun()
            for i in range(len(self.deps)):
                if i < curdep:
                    self.deps[i].startDiscreteRun()
        for i, lb in enumerate(self.lbs):
            if lb['do_linking']:
                assert self.headnode == lb['lb'].dep_idx
                # builtins.loglinecount += 1
                # logger.error("Checkpoint # %d" % builtins.loglinecount)
                # if self.G.number_of_nodes() == 16:
                #     logger.error("J")
                lb['lb'].connect_cycles(self.DG.nodes[self.headnode]['graph'], lb['global_graph'], self.G)

        # if builtins.loglinecount == 30 and self.ord == 3 and self.idx==0:
        #     self.masterunit.masterunit.double_check(hard=True)

        self.DG.nodes[self.headnode]['ikprob'].updatexyz(self.G)

        # if builtins.loglinecount == 30 and self.ord == 3 and self.idx==0:
        #     self.masterunit.masterunit.double_check(hard=True)

        if getbool("DoValidation", "IK_CyclicPart", self.config):
            # xyzlines = [str(self.G.number_of_nodes()), ""]
            # for node in list(self.G.nodes()):
            #     xyzlines.append("%3s%10.4f%10.4f%10.4f%3s" % ("C", self.G.nodes[node]['xyz'][0],
            #                                                   self.G.nodes[node]['xyz'][1],
            #                                                   self.G.nodes[node]['xyz'][2], ""))
            # wfile = open("indepframe.xyz", "w")
            # wfile.write("\n".join(xyzlines))
            # wfile.close()

            good = self.validator.validate(self.G)
            logger.info("Validator returned " + repr(good))
            if not good:
                raise Exception("Failed to obtain IK solution")

        for cnode_idx, cnode in enumerate(self.G.nodes):
            logger.debug("%d node before: %s" % (cnode, repr(self.G.nodes[cnode]['xyz'])))

            myframe = None
            for frame in self.G.nodes[cnode]['frames']:
                if self.headnode in frame['users']:
                    myframe = frame
                    break
            assert myframe is not None

            at1 = self.G.nodes[cnode]['xyz']
            at0 = self.G.nodes[myframe['innerframe'][0]]['xyz']
            at2 = self.G.nodes[myframe['innerframe'][1]]['xyz']
            self.extG.nodes[cnode]['xyz'] = deepcopy(at1) # NOT NECESSARILY COPY

            xv = at0 - at1
            yv = at2 - at1
            zv = IK_Math.gs_rand(xv, yv)

            for i in range(len(myframe['outerNB'])):
                nbatom = myframe['outerNB'][i]
                if nbatom not in self.G.nodes():
                    logger.debug("Setting position of %d atom" % myframe['outerNB'][i])
                    localcoord = myframe['localframeXYZ'][i]
                    self.extG.nodes[nbatom]['xyz'] = at1 + \
                                                     localcoord[0] * xv + \
                                                     localcoord[1] * yv + \
                                                     localcoord[2] * zv
                    myframe['outerXYZ'][i] = self.extG.nodes[nbatom]['xyz']
            logger.debug("%d node after: %s" % (cnode, repr(self.G.nodes[cnode]['xyz'])))

        logger.info("Exiting GU order=%d; idx=%d" % (self.ord, self.idx))
        if getbool("RunGeomCheck", "IK_GeomUnit", self.config):
            if geomcheckTree_graph(self.extG, getfloat("MinDistance", "IK_GeomUnit", self.config)):
                logger.info("Geomcheck passed")
                self.FullPS.cause = cause.success
                if self.scounter.hasddof:
                    self.good_ddof = self.scounter.ddof.getValue()
                    logger.info("Setting gddof=%d order=%d; idx=%d" % (self.good_ddof, self.ord, self.idx))
                return IK_ParameterSet()
            else:
                logger.info("Geomcheck failed on node %d (order=%d; idx=%d)" % (self.headnode, self.ord, self.idx))
                self.FullPS.cause = cause.geomoverlap
                if increase_discrete and not incremented:
                    if not self.scounter.hasddof:
                        raise Exception("Can't fix that overlap!")
                    logger.info("Reading gddof=%d order=%d; idx=%d" % (self.good_ddof, self.ord, self.idx))
                    self.scounter.ddof.setValue(self.good_ddof)
                    self.scounter.value = self.good_ddof
                    return self.DG.nodes[self.headnode]['ikprob'].getPS(errorcause=cause.geomoverlap_ddof)
                logger.info("Excel: line 257 resetDDOF")
                self.scounter.resetDDOF()
                return self.DG.nodes[self.headnode]['ikprob'].getPS(errorcause=cause.geomoverlap)
        else:
            self.FullPS.cause = cause.success
            return IK_ParameterSet()

    def transmit_graphs(self, desired_node, source=None, trymaster=False):
        if desired_node == self.headnode:
            return (self.DG.nodes[self.headnode]['graph'], self.G)
        else:
            if desired_node in self.dep_indicies:
                return self.deps[self.dep_indicies.index(desired_node)].transmit_graphs(desired_node)

            for dep in self.deps:
                if dep.headnode == source:
                    continue
                res = dep.transmit_graphs(desired_node)
                if res is not None:
                    return res

            if trymaster:
                if self.masterunit is not None:
                    return self.masterunit.transmit_graphs(desired_node, source=self.headnode, trymaster=True)
                else:
                    raise Exception("Can not find the desired node")
            else:
                return None

    def double_check(self, hard=False):
        for dep in self.deps:
            dep.double_check()
        if len(self.deps) == 0 or hard:
            xyzlines = [str(self.G.number_of_nodes()), ""]
            for node in list(self.G.nodes()):
                xyzlines.append("%3s%10.4f%10.4f%10.4f%3s" % ("C", self.G.nodes[node]['xyz'][0],
                                                              self.G.nodes[node]['xyz'][1],
                                                              self.G.nodes[node]['xyz'][2], ""))
            wfile = open("indepframe.xyz", "w")
            wfile.write("\n".join(xyzlines))
            wfile.close()
            self.validator.validate(self.G, weak=True)

    def perturb_geometry(self):
        for dep in self.deps:
            dep.perturb_geometry()
        self.DG.nodes[self.headnode]['ikprob'].perturb_geometry()
