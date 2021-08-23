from copy import copy
import numpy as np
from numpy.linalg import norm
from enum import Enum

from .math import IK_Math
from .logging import createLogger

logger = createLogger("DataStructures")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class itemtype(Enum):
    BOND = 1
    ATOM = 2
    UNDEFINED = 3

class Cngroup:
    def __init__(self, size, itype=itemtype.UNDEFINED):
        self.ord = size
        self.itype = itype
        self.items = []
        for i in range(self.ord):
            newinst = GroupItem(self)
            self.items.append(newinst)
            self.items[i]['num'] = i

    def atoms_find_proper_atoms(self):
        res = []
        for atom in self.items:
            if atom['bonds'][0]['free'] and atom['bonds'][1]['free']:
                res.append(atom)
        return res

    def atoms_dist(self, a, b):
        dist = abs(a['num'] - b['num']) - 1
        return min(self.ord - dist - 2, dist)

    def bond_getstart(self):
        for item in self.items:
            if not item['free']:
                return item + 1
        return self.items[0]

    def atoms_getstep(self, pair, axatoms):
        curatom = pair[0]
        otherat = copy(axatoms)
        otherat.remove(pair[0])
        otherat.remove(pair[1])
        otherat = otherat[0]
        while not curatom == pair[1]:
            curatom += 1
            if curatom == otherat:
                return -1
        return 1

    def __len__(self):
        return self.ord

    def __getitem__(self, key):
        return self.items[key % self.ord]

    def perturb_angles(self, sd):
        for item in self.items:
            item.vangle = np.random.normal(item.vangle_original, sd * deg2rad)

    def perturb_bonds(self, sd):
        for item in self.items:
            item.length = np.random.normal(item.length_original, sd)

class GroupItem:
    def __init__(self, group):
        self.group = group
        if group.itype == itemtype.ATOM:
            self.vangle = None
        elif group.itype == itemtype.BOND:
            self.length = None

    def bond_getsideatoms(self):
        return [(self['atoms'][0]-1)['G_idx'], (self['atoms'][1]+1)['G_idx']]

    def bondbetween(self, other):
        if self.group.itype == itemtype.ATOM and other.group.itype == itemtype.ATOM:
            if self + 1 == other:
                if other in self['bonds'][1]['atoms']:
                    return self['bonds'][1]
            elif other + 1 == self:
                if other in self['bonds'][0]['atoms']:
                    return self['bonds'][0]
            else:
                raise Exception('Atoms must be neighbours!')
        else:
            raise Exception('Only for atoms!')
        raise Exception('Something bad happened')

    def bond_gettorsion(self, attr="xyz"):
        at_p = self['atoms'][0]
        at_n = self['atoms'][1]
        side_p = (at_p - 1)[attr]
        side_n = (at_n + 1)[attr]
        at_p = at_p[attr]
        at_n = at_n[attr]
        return IK_Math.gettorsion([side_p, at_p, at_n, side_n])

    def atom_getangle(self, forcecalc=False, attr="xyz"):
        if self.vangle == None or forcecalc:
            prevvec = (self - 1)[attr]
            curvec = self[attr]
            nextvec = (self + 1)[attr]
            dir1 = prevvec - curvec
            dir2 = nextvec - curvec
            dir1 /= norm(dir1)
            dir2 /= norm(dir2)
            angle = np.arccos(np.dot(dir1, dir2))
            if self.vangle == None:
                self.vangle = angle
                # TODO if getbool("AllowAnglePerturbation", "IK_TLCSolver", self.config):
                self.vangle_original = angle
                logger.debug("Vangle %d(=%d) is set" % (self.group.items.index(self), self['G_idx']))
            return angle
        return self.vangle

    def bond_getcommonatom(self, other):
        if self['atoms'][0] in other['atoms']:
            return self['atoms'][0]
        elif self['atoms'][1] in other['atoms']:
            return self['atoms'][1]

    def bond_getlength(self, forcecalc=False, attr="xyz"):
        if self.length == None or forcecalc:
            result = norm(self['atoms'][0][attr] - self['atoms'][1][attr])
            if self.length == None:
                self.length = result
                # TODO if getbool("AllowBondPerturbation", "IK_TLCSolver", self.config):
                self.length_original = result
            return result
        return self.length

    def __add__(self, other):
        if isinstance(other, int):
            idx = other
        elif isinstance(other, GroupItem):
            idx = other['num']
        return self.group[idx + self['num']]

    def __sub__(self, other):
        if isinstance(other, int):
            idx = other
        elif isinstance(other, GroupItem):
            idx = other['num']
        return self.group[self['num'] - idx]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __eq__(self, other):
        return id(self.group) == id(other.group) and self['num'] == other['num']

    def __repr__(self):
        if self.group.itype == itemtype.BOND:
            return "Bond #%d : %d(=%d) - %d(=%d)" % (self['num'], self['atoms'][0]['num'], self['atoms'][0]['G_idx'], \
                                                     self['atoms'][1]['num'], self['atoms'][1]['G_idx'])
        elif self.group.itype == itemtype.ATOM:
            return "Atom %d(=%d)" % (self['num'], self['G_idx'])
        else:
            return "Object of undefined type #%d" % (self['num'])
