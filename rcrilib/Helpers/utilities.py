from copy import deepcopy
from enum import Enum

from .datastructures import GroupItem
from .logging import createLogger

logger = createLogger("Utilities")

class ikdof(Enum):
    CONTINUOUS = 1
    DISCRETE = 2
    FREE = 3
    DEPENDENT = 4
    FIXED = 5

class target_class(Enum):
    MOLECULE = 1
    CYCLICPART = 2
    PROBLEM = 3
    SOLVER = 4
    RIGIDFRAG = 5
    FOURFIVE = 6

class cause(Enum):
    unknown = 0
    zerosolutions = 1
    geomoverlap = 2
    zerosolutions_ddof = 3
    geomoverlap_ddof = 4
    success = 5

class IK_ParameterSet:  # NEED TO OVERLOAD + TO ADD NEW PARAMETERS OF THE CYCLE
    def __init__(self, p = None):
        self.params = []
        self.cause = cause.unknown
        self.success = True
        if p is not None:
            self += p

    def __repr__(self):
        nfixed = 0
        ndep = 0
        nfree = 0
        nnoncyc = 0
        for i in self.params:
            if i.isFixed():
                nfixed += 1
            elif i.isDependent():
                ndep += 1
            elif i.isFree():
                nfree += 1
        return "Parameter set. Total: %d free parameters, %d dependent, %d fixed.\n List of parameters: %s" % (
            nfree, ndep, nfixed, repr(self.params))

    def __iadd__(self, other):
        if isinstance(other, IK_Parameter):
            self.params.append(other)
        elif isinstance(other, IK_ParameterSet):
            for par in other.params:
                self.params.append(par)
        return self

    def __len__(self):
        return len(self.params)

    def getPS(self, excludeFixed=False, excludeFree=False, excludeDep=False, excludeDiscrete=False, errorcause=cause.unknown):
        if errorcause != cause.unknown:
            excludeFixed = True
            excludeDep = True
        newPS = IK_ParameterSet()
        for item in self.params:
            if not (excludeFixed and item.isFixed() or excludeDep and item.isDependent() or
                    excludeFree and item.isFree() or excludeDiscrete and item.dim == ikdof.DISCRETE):
                newPS.params.append(item)
        newPS.cause = errorcause
        if self.cause != cause.unknown and errorcause == cause.unknown:
            newPS.cause = self.cause
        if len(newPS.params) > 0 or not self.success:
            newPS.success = False
        if errorcause != cause.success and errorcause != cause.unknown:
            newPS.success = False
        return newPS

    def byTarget(self, tc, errorcause=cause.unknown):
        newPS = IK_ParameterSet()
        for item in self.params:
            if item.tc == tc:
                newPS.params.append(item)
        newPS.cause = errorcause
        if len(newPS) > 0 or errorcause != cause.success:
            newPS.success = False
        return newPS

    def __iter__(self):
        return iter(self.params)

    def __getitem__(self, key):
        return self.params[key]

    def isEmpty(self):
        return len(self.params) == 0

    def isChanged(self):
        if len(self.params) == 0:
            return True
        for param in self.params:
            if param.changed:
                return True
        return False

    def freeze(self):
        for param in self.params:
            param.changed = False

    def hasDDOF(self):
        for param in self.params:
            if param.isDiscrete():
                return True
        return False

    def getMyDiscreteParameter(self, tc):
        for param in self.params:
            if param.tc == tc and param.isDiscrete():
                return param

def dofstr(mytype):
    if mytype == ikdof.FREE:
        return "Free"
    elif mytype == ikdof.DEPENDENT:
        return "Dependent"
    elif mytype == ikdof.FIXED:
        return "Fixed"
    elif mytype == ikdof.NONCYCLIC:
        return "Non-cyclic"

class IK_Parameter:
    def __init__(self, dim, dep, sideatoms, tc, bond=None, atoms=None):
        self.dim = dim
        self.dep = dep
        self.tc = tc
        self.sides = deepcopy(sideatoms)

        if dim == ikdof.CONTINUOUS:
            if tc == target_class.MOLECULE:
                self.frags = []

            if isinstance(bond,GroupItem):
                self.bond = bond
                self.atoms = [bond['atoms'][0]['G_idx'],bond['atoms'][1]['G_idx']]
            else:
                self.atoms = deepcopy(atoms)
        self.value = None
        if dim == ikdof.DISCRETE:
            self.maxValue = None
        self.changed = True

    def __repr__(self):
        if self.dim == ikdof.CONTINUOUS:
            return "%s parameter of %d-%d bond" % (dofstr(self.dep), self.atoms[0], self.atoms[1])
        elif self.dim == ikdof.DISCRETE:
            return "Discrete degree of freedom"

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def isDiscrete(self):
        return self.dim == ikdof.DISCRETE

    def isContinuous(self):
        return self.dim == ikdof.CONTINUOUS

    def isFree(self):
        return self.dep == ikdof.FREE

    def isDependent(self):
        return self.dep == ikdof.DEPENDENT

    def isFixed(self):
        return self.dep == ikdof.FIXED

    def setValue(self, newval):
        if self.value != newval:
            self.value = newval
            self.changed = True

    def setCounter(self, sc):
        if self.dim == ikdof.DISCRETE:
            self.solutionCounter = sc
        else:
            raise Exception("Fix it")

# TODO Try to make these three @chached
def getfloat(myoption, myclass, config):
    if config.has_option(myclass, myoption):
        return config[myclass].getfloat(myoption)
    else:
        return config["IK_All"].getfloat(myoption)

def getint(myoption, myclass, config):
    if config.has_option(myclass, myoption):
        return config[myclass].getint(myoption)
    else:
        return config["IK_All"].getint(myoption)

def getbool(myoption, myclass, config):
    if config.has_option(myclass, myoption):
        return config[myclass].getboolean(myoption)
    else:
        return config["IK_All"].getboolean(myoption)
