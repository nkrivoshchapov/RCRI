from .logging import createLogger
import math, multiprocessing

logger = createLogger("SystematicSearch")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

class JointDriver:
    def __init__(self,cyc,noncyc):
        self.cycdriver = cyc
        self.noncycdriver = noncyc
        logger.info("Constructed joint syssearch driver.")

    def isDone(self):
        return self.cycdriver.isDone() and self.noncycdriver.isDone()

    def switchToNext(self):
        if self.noncycdriver.isDone():
            self.noncycdriver.reset()
            self.cycdriver.switchToNext()
        else:
            self.noncycdriver.switchToNext()

    def setValues(self):
        self.cycdriver.setValues()
        self.noncycdriver.setValues()

    def reset(self):
        self.cycdriver.reset()
        self.noncycdriver.reset()

class SystematicDriver:
    def __init__(self, ps, stepsize, startvalue):
        self.PS = ps
        logger.info("Creating syssearch driver. PS = " + repr(self.PS))

        self.cursteps = [0] * len(self.PS.params)
        self.stepsize = stepsize
        self.startvalue = startvalue
        self.maxvalue = math.ceil(360/stepsize) - 1

    def isDone(self):
        for item in self.cursteps:
            if item != self.maxvalue:
                return False
        return True

    def setValues(self):
        for i, item in enumerate(self.cursteps):
            self.PS[i].setValue((item*self.stepsize + self.startvalue)*deg2rad)

    def switchToNext(self):
        increased = False
        for i, item in enumerate(self.cursteps):
            if item == self.maxvalue:
                self.cursteps[i] = 0
            else:
                self.cursteps[i] += 1
                increased = True
                break
        if not increased:
            raise Exception("Cannot make another step. Systematic search is finished")

    def reset(self):
        # if not self.isDone():
        #     raise Exception("Search is not finished!")
        self.cursteps = [0] * len(self.PS.params)

class SystematicDriver_parallel(object):
    def __init__(self, ps, stepsize, startvalue):
        super(SystematicDriver_parallel, self).__init__()
        self.mgr = multiprocessing.Manager()
        self.lock = multiprocessing.Lock()

        logger.info("Creating parallelized syssearch driver. PS = " + repr(ps))
        self.cursteps = self.mgr.list([0] * len(ps.params))
        self.stepsize = stepsize
        self.startvalue = startvalue
        self.maxvalue = math.ceil(360 / stepsize) - 1

    def switchToNext(self):
        with self.lock:
            increased = False
            for i, item in enumerate(self.cursteps):
                if item == self.maxvalue:
                    self.cursteps[i] = 0
                else:
                    self.cursteps[i] += 1
                    increased = True
                    break
            if not increased:
                raise Exception("Cannot make another step. Systematic search is finished")

    def isDone(self):
        with self.lock:
            for item in self.cursteps:
                if item != self.maxvalue:
                    return False
            return True

    def setValues(self, ps):
        with self.lock:
            for i, item in enumerate(self.cursteps):
                ps[i].setValue((item*self.stepsize + self.startvalue)*deg2rad)

class Counter_parallel(object):
    def __init__(self, initval=0):
        self.val = multiprocessing.Value('i', initval)
        self.lock = multiprocessing.Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value
