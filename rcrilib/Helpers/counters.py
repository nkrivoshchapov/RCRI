from copy import copy
import builtins

from .logging import createLogger
from .utilities import cause, ikdof, target_class as tc

logger = createLogger("IK_Counters")

class IK_SolutionCounter:
    def __init__(self, ps):
        self.ddof = ps.getMyDiscreteParameter(tc.SOLVER)
        if self.ddof is not None:
            self.hasddof = True
            ps.getMyDiscreteParameter(tc.SOLVER).setCounter(self)  # TODO What if no DDOF in ps?
        else:
            self.hasddof = False

        self.ps = ps.getPS(excludeDiscrete=True)

        self.value = None
        self.dep = []
        for param in self.ps:
            if param.dep == ikdof.DEPENDENT:
                self.dep.append(copy(param.atoms))
        self.cont_dofs = []
        for i in self.ps.params:
            self.cont_dofs.append(None)

    def recordState(self):
        if self.hasddof:
            self.value = self.ddof.value
            self.maxValue = self.ddof.maxValue
            for i in range(len(self.cont_dofs)):
                self.cont_dofs[i] = self.ps.params[i].value
            self.checkState()

    def checkState(self):
        if self.hasddof:
            if self.ddof.value is not None and self.value != self.ddof.value:
                raise Exception("Unexpected change of DDOF value %s vs %s" % (repr(self.value), repr(self.ddof.value)))

            if self.maxValue != self.ddof.maxValue:
                raise Exception("Unexpected change of DDOF maxValue")

            if self.ddof.value is not None:
                for i in range(len(self.cont_dofs)):
                    if self.cont_dofs[i] != self.ps.params[i].value:
                        raise Exception("Unexpected change of continuous DOF value")

    def checkPS(self, ps): # TODO Remove this functions?
        if self.hasddof:
            if hasattr(self,"maxValue"):
                self.checkState()
            if self.ddof.value is not None and ps.cause != cause.geomoverlap_ddof and ps.cause != cause.zerosolutions_ddof:
                bad = False
                for param in ps.params:
                    if param.dim == ikdof.DISCRETE:
                        continue
                    if param in self.ps.params or param.atoms in self.dep:
                        logger.error("Dont change " + repr(param))
                        bad = True
                if bad:
                    raise Exception("Change of parameter in PS is going to affect the number of IK solutions!")

    def isDepChanged(self, ps):
        if self.hasddof:
            if hasattr(self,"maxValue"):
                self.checkState()
            if self.ddof.value is not None and ps.cause != cause.geomoverlap_ddof and ps.cause != cause.zerosolutions_ddof:
                for param in ps.params:
                    if param.dim == ikdof.DISCRETE:
                        continue
                    if param.atoms in self.dep:
                        logger.error("Dont change " + repr(param))
                        return True
            return False
        else:
            return False

    def startDiscreteRun(self):
        if self.hasddof:
            self.startDiscr = self.ddof.value
            self.maxDiscr = self.ddof.maxValue
            self.discValues = list(range(self.maxDiscr))
            self.discValues.remove(self.startDiscr)

    def done_discr(self):
        if self.hasddof:
            return len(self.discValues) == 0
        else:
            return True

    def nextValue(self):
        if self.hasddof:
            logger.info("Excel: Increased DDOF")
            logger.info("Short: Increased DDOF")
            self.value = self.discValues.pop()
            self.ddof.setValue(self.value)

    def resetDDOF(self):
        if self.hasddof:
            self.discValues = []
            self.ddof.setValue(None)

    def logstate(self, nocount = False):
        if not self.hasddof:
            return "%s;NoDDOF" % (self.name)

        if hasattr(self, "discValues"):
            if self.ddof.value in self.discValues:
                raise Exception("Repeated DDOFs!!!")
        if nocount:
            if hasattr(self, "discValues"):
                return "%s;%s;%s" % (self.name, repr(self.discValues), repr(self.ddof.value))
            else:
                return "%s;NAN;%s" % (self.name, repr(self.ddof.value))
        else:
            builtins.scounter_linecount += 1
            if hasattr(self,"discValues"):
                return "%s;%s;%s;line #%d" % (self.name, repr(self.discValues), repr(self.ddof.value), builtins.scounter_linecount)
            else:
                return "%s;NAN;%s;line #%d" % (self.name, repr(self.ddof.value), builtins.scounter_linecount)

    def backup_ddof(self):
        if self.hasddof:
            self.startDiscr_backup = self.startDiscr
            self.maxDiscr_backup = self.maxDiscr
            self.discValues_backup = copy(self.discValues)
            self.value_backup = self.value

    def restore_ddof(self):
        if self.hasddof:
            self.startDiscr = self.startDiscr_backup
            self.maxDiscr = self.maxDiscr_backup
            self.discValues = copy(self.discValues_backup)
            self.value = self.value_backup
            self.ddof.maxValue = self.maxDiscr_backup
            self.ddof.setValue(self.value_backup)

class IK_CompoundCounter:
    def __init__(self):
        self.counters = []

    def __iadd__(self, other):
        if isinstance(other, IK_SolutionCounter):
            self.counters.append(other)
        elif isinstance(other, IK_CompoundCounter):
            for item in other.counters:
                self.counters.append(item)
        else:
            raise Exception("Fix it pls")
        return self

    def recordState(self):
        for item in self.counters:
            item.recordState()

    def checkState(self):
        for item in self.counters:
            item.checkState()

    def checkPS(self, ps):
        for item in self.counters:
            item.checkPS(ps)

    def startDiscreteRun(self):
        for item in self.counters:
            item.startDiscreteRun()

    def done_discr(self):
        for item in self.counters:
            if not item.done_discr():
                return False
        return True

    def resetDDOF(self):
        for item in self.counters:
            item.resetDDOF()

    def logstate(self, nocount = False):
        parts = []
        for item in self.counters:
            parts.append(item.logstate(nocount))
        return " // ".join(parts)

    def backup_ddof(self):
        for item in self.counters:
            item.backup_ddof()

    def restore_ddof(self):
        for item in self.counters:
            item.restore_ddof()