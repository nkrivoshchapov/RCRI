import configparser, builtins, os, traceback, sys, glob

def single_threaded(configfile, testing = False):
    config = configparser.ConfigParser()
    config.read(configfile)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    builtins.scounter_linecount = 0

    from rcrilib import singthr_mc, Timings
    from rcrilib.Helpers import createLogger

    timer = Timings()
    logger = createLogger("MainScript")
    logger.error("Starting IK-MC conformer generation with bruteforce of DDOFs")

    if testing:
        todolist = glob.glob(config["Testing"].get("TestMolecules"))
        start = 0
        skipfile = config["Testing"].get("SkipTill")
        if len(skipfile) > 0 and skipfile in todolist:
            start = todolist.index(skipfile)
        logger.error("Skipped %d files" % start)
        for i in range(start+1, len(todolist)):
            singthr_mc(todolist[i], config, logger, timer, testing)
            timer.summarize(todolist[i])
            timer.reset()
    else:
        molfile = config["MCDriver"].get("InputFile")
        if not os.path.isfile(molfile):
            raise Exception("Input file \"%s\" is not found" % molfile)
        singthr_mc(molfile, config, logger, timer, testing)

def parallel(configfile, numthreads, testing = False):
    config = configparser.ConfigParser()
    config.read(configfile)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    builtins.scounter_linecount = 0

    from rcrilib import run_parallel_mc, Timings_parallel
    from rcrilib.Helpers import createLogger

    timer = Timings_parallel()
    logger = createLogger("MainScript")
    logger.warning("Starting IK-MC conformer generation with bruteforce of DDOFs")

    if testing:
        todolist = glob.glob(config["Testing"].get("TestMolecules"))
        start = 0
        skipfile = config["Testing"].get("SkipTill")
        if len(skipfile) > 0 and skipfile in todolist:
            start = todolist.index(skipfile)
        logger.error("Skipped " + repr(start))
        for i in range(start+1, len(todolist)):
            run_parallel_mc(todolist[i], config, configfile, numthreads, logger, timer, testing)
            timer.summarize(todolist[i])
    else:
        molfile = config["MCDriver"].get("InputFile")
        if not os.path.isfile(molfile):
            raise Exception("Input file \"%s\" is not found" % molfile)
        run_parallel_mc(molfile, config, configfile, numthreads, logger, timer, testing)

if __name__ == '__main__':
    if "-c" in sys.argv:
        configfile = sys.argv[sys.argv.index("-c")+1]
    else:
        configfile = "default.ini"
    print("Using config-file " + configfile)

    if "-T" in sys.argv:
        numthreads = int(sys.argv[sys.argv.index("-T")+1])
    else:
        numthreads = 1

    testing = False
    if "--test" in sys.argv:
        testing = True
        print("--test -> Run testing")

    if numthreads == 1:
        print("Running single-threaded version")
        single_threaded(configfile, testing=testing)
    else:
        print("Running in %d threads" % numthreads)
        parallel(configfile, numthreads, testing=testing)

    # parallel("default.ini", 4, testing=True)
