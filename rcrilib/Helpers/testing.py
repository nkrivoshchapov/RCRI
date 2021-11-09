import time, multiprocessing

class Timings:
    def __init__(self):
        self.goodtimes = []
        self.badtimes = []
        self.inittimes = []

    def reset(self, man=None):
        self.goodtimes = []
        self.badtimes = []
        self.inittimes = []

    def check_start(self):
        self.start = time.perf_counter()

    def stop_init(self):
        # stoptime = time.perf_counter()
        self.inittimes.append(time.perf_counter())
        # self.start = stoptime

    def stop_good(self):
        # stoptime = time.perf_counter()
        self.goodtimes.append(time.perf_counter())
        # self.start = stoptime

    def stop_bad(self):
        # stoptime = time.perf_counter()
        self.badtimes.append(time.perf_counter())
        # self.start = stoptime

    def summarize(self, molfile):
        for i in range(len(self.goodtimes)):
            self.goodtimes[i] -= self.start
        for i in reversed(range(1, len(self.goodtimes))):
            self.goodtimes[i] -= self.goodtimes[i-1]

        for i in range(len(self.badtimes)):
            self.badtimes[i] -= self.start
        for i in reversed(range(1, len(self.badtimes))):
            self.badtimes[i] -= self.badtimes[i-1]

        for i in range(len(self.inittimes)):
            self.inittimes[i] -= self.start
        for i in reversed(range(1, len(self.inittimes))):
            self.inittimes[i] -= self.inittimes[i-1]

        print("-" * 40)
        print("| Timings for %s" % molfile)
        print("|%34s : %f s" % (
            "Average time for initialization", sum(self.inittimes) / len(self.inittimes)))
        if len(self.goodtimes) > 0:
            print("|%34s : %d (average time: %f s)" % (
                "Number of successful iterations", len(self.goodtimes), sum(self.goodtimes) / len(self.goodtimes)))
        else:
            print("|%34s : %d" % (
                "Number of successful iterations", len(self.goodtimes)))
        if len(self.badtimes) > 0:
            print("|%34s : %d (average time: %f s)" % (
                "Number of unsuccessful iterations", len(self.badtimes), sum(self.badtimes) / len(self.badtimes)))
        else:
            print("|%34s : %d" % (
                "Number of unsuccessful iterations", len(self.badtimes)))
        print("-" * 40)

class Timings_parallel(Timings):
    def reset(self, man=None):
        self.goodtimes = man.list()
        self.badtimes = man.list()
        self.inittimes = man.list()

    def format_lists(self):
        self.goodtimes = list(self.goodtimes)
        self.badtimes = list(self.badtimes)
        self.inittimes = list(self.inittimes)
