from .logging import createLogger

logger = createLogger("sharedValue")

class sharedValue:
    def __init__(self, descr):
        self.descr = descr
        self.value = None
        self.changed = True
        logger.info("Creating shared value: " + self.descr)

    def __repr__(self):
        if isinstance(self.value, sharedValue):
            return "Dependent shared value - " + self.descr
        else:
            return "Shared value - " + self.descr

    def getValue(self):
        if isinstance(self.value, sharedValue):
            return self.value.getValue()
        else:
            return self.value

    def setValue(self, v, trick=False):
        if isinstance(v, sharedValue):
            raise Exception("Shouldn't assign a shared value")
        if isinstance(self.value, sharedValue):
            raise Exception("Shouldn't assign value to dependent shared value")
        self.value = v
        if not trick:
            self.changed = True

    def freeze(self):
        self.changed = False

    def isChanged(self):
        return self.changed

    def __float__(self):
        return self.getValue()
