import time

class Testing:
    def __init__(self):
        self.goodtimes = []
        self.badtimes = []
        self.inittimes = []

    def reset(self):
        self.goodtimes = []
        self.badtimes = []
        self.inittimes = []

    def check_start(self):
        self.start = time.perf_counter()

    def stop_init(self):
        stoptime = time.perf_counter()
        self.inittimes.append(time.perf_counter() - self.start)
        self.start = stoptime

    def stop_good(self):
        stoptime = time.perf_counter()
        self.goodtimes.append(time.perf_counter() - self.start)
        self.start = stoptime

    def stop_bad(self):
        stoptime = time.perf_counter()
        self.badtimes.append(time.perf_counter() - self.start)
        self.start = stoptime

    def summarize(self, molfile):
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
