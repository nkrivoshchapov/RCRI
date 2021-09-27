import configparser, builtins, random, os, sys, multiprocessing, glob
from copy import deepcopy

def writeGeom(config, mol, logger, conf_count):
    if ".mol" in config["MCDriver"].get("OutputFile"):
        mol.writeToMol(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
        logger.warning("Writing file " + config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    elif ".xyz" in config["MCDriver"].get("OutputFile"):
        mol.writeToXyz(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
        logger.warning("Writing file " + config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    logger.info("Conformer was generated successfully.")

def termination_crit(count, maxnumber):
    if maxnumber == -1:
        return False
    else:
        return count >= maxnumber

def worker(s, pool, conf_count, molfile, configname, testing):
    procname = multiprocessing.current_process().name
    config = configparser.ConfigParser()
    config.read(configname)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    builtins.scounter_linecount = 0

    from rcrilib import IK_Molecule, SystematicDriver
    from rcrilib.Helpers import createLogger, log_versions, cause
    do_syssearch_noncyclic = config["SystematicSearch"].getboolean("EnableNonCyclic")
    do_syssearch_cyclic = config["SystematicSearch"].getboolean("EnableCyclic")

    logger = createLogger("WorkerScript")
    logger.warning("Starting IK-MC conformer generation with bruteforce of DDOFs")
    if not os.path.isfile(molfile):
        raise Exception("Input file \"%s\" is not found" % molfile)

    if testing:
        maxconf = config["Testing"].getint("NumberOfConfs")
    else:
        maxconf = config["MCDriver"].getint("NumberOfConfs")
    log_versions(logger)
    logger.warning("Initializing from file \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    if do_syssearch_cyclic:
        sysdriver_cyc = pool
    if do_syssearch_noncyclic:
        sysdriver_noncyc = SystematicDriver(mol.PS.getNonCyclicDOFs(),
                                            config["SystematicSearch"].getfloat("NonCyclicStep"),
                                            config["SystematicSearch"].getfloat("NonCyclicStart"))

    while not termination_crit(conf_count.value(), maxconf):
        logger.warning("Generated %d conformations. Goal = %d" % (conf_count.value(), maxconf))
        ps = mol.getPS()
        ps.success = False
        mol.discr_cp = -1
        it = 0
        logger.warning("Generating the initial conformation")
        while it != config["MCDriver"].getint("MaxTries") and not ps.success:
            for item in ps:
                if item.isContinuous():
                    item.setValue(random.uniform(-3.141592, 3.141592))
                if item.isDiscrete():
                    item.setValue(None)
            if it == 0 and do_syssearch_cyclic:
                with s:
                    if not sysdriver_cyc.isDone():
                        logger.warning("Making a step in cyclic DOFs")
                        sysdriver_cyc.setValues(mol.PS.getCyclicDOFs())
                        sysdriver_cyc.switchToNext()
            ps = mol.applyPS()
            it += 1
            if ps.success:
                if not testing:
                    writeGeom(config, mol, logger, conf_count.value())
                conf_count.increment()
                break
            elif ps.isEmpty():
                it = config["MCDriver"].getint("MaxTries")
                break
        if it == config["MCDriver"].getint("MaxTries"):
            mol.perturb_geometry()
            logger.warning("Couldn't generate the initial conformation. Trying again.")
            continue
        elif termination_crit(conf_count.value(), maxconf):
            break

        backup_ps = deepcopy(mol.getPS())
        if do_syssearch_noncyclic:
            logger.warning("Starting non-cyclic search")
            while not sysdriver_noncyc.isDone() and \
                                not termination_crit(conf_count.value(), maxconf):
                sysdriver_noncyc.setValues()
                sysdriver_noncyc.switchToNext()
                ps = mol.applyPS()
                if ps.success:
                    if not testing:
                        writeGeom(config, mol, logger, conf_count.value())
                    conf_count.increment()
            sysdriver_noncyc.reset()
            ps = mol.getPS()
            for i in range(len(ps)):
                if ps[i].isContinuous():
                    if ps[i].value != backup_ps[i].value:
                        ps[i].setValue(backup_ps[i].value)

        mol.startDiscreteRun()
        trycount = 0
        prev = ""
        logger.warning("Starting bruteforce of DDOF combinations")
        while not mol.done_all_discr() and \
                                not termination_crit(conf_count.value(), maxconf):
            logger.warning("New DDOF combination")
            check = mol.ccounter.logstate(nocount=True)
            ps = mol.applyPS(increase_discrete=True)
            mol.ccounter.checkState()
            if prev == check:
                raise Exception("DDOF state repetition!!!")
            prev = check
            trycount += 1
            it = 0
            while it != config["MCDriver"].getint("MaxTries") and not ps.success and \
                    ps.cause != cause.geomoverlap_ddof and ps.cause != cause.zerosolutions_ddof:
                for item in ps:
                    if item.isContinuous():
                        logger.info("Setting random torsion on " + repr(item.atoms))
                        item.setValue(random.uniform(-3.141592, 3.141592))
                ps = mol.applyPS()
                it += 1
                if not ps.success and ps.isEmpty():
                    break
            if ps.success:
                if not testing:
                    writeGeom(config, mol, logger, conf_count.value())
                conf_count.increment()
                if termination_crit(conf_count.value(), maxconf):
                    break

            if do_syssearch_noncyclic:
                logger.warning("Starting non-cyclic search")
                while not sysdriver_noncyc.isDone() and \
                                not termination_crit(conf_count.value(), maxconf):
                    sysdriver_noncyc.setValues()
                    sysdriver_noncyc.switchToNext()
                    ps = mol.applyPS()
                    if ps.success:
                        if not testing:
                            writeGeom(config, mol, logger, conf_count.value())
                        conf_count.increment()
                sysdriver_noncyc.reset()

            ps = mol.getPS()
            for i in range(len(ps)):
                if ps[i].isContinuous():
                    if ps[i].value != backup_ps[i].value:
                        ps[i].setValue(backup_ps[i].value)
            mol.ccounter.recordState()
        logger.warning("DDOF bruteforce is finished")
        mol.perturb_geometry()
        if do_syssearch_cyclic and config["SystematicSearch"].getboolean("TerminateWhenDone"):
            with s:
                if sysdriver_cyc.isDone():
                    logger.warning("Systematic search is finished. Terminating.")
                    break

def main(configfile, numthreads, testing = False):
    config = configparser.ConfigParser()
    config.read(configfile)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    builtins.scounter_linecount = 0
    from rcrilib import IK_Molecule, SystematicDriver_parallel, Counter_parallel, Testing
    from rcrilib.Helpers import createLogger
    logger = createLogger("MainScript")
    molfile = config["MCDriver"].get("InputFile")

    def run_parallel_mc(molfile):
        if not os.path.isfile(molfile):
            raise Exception("Input file \"%s\" is not found" % molfile)

        mol = IK_Molecule(molfile, config)
        mol.prepare_ik()
        pool = SystematicDriver_parallel(mol.PS.getCyclicDOFs(),
                                         config["SystematicSearch"].getfloat("CyclicStep"),
                                         config["SystematicSearch"].getfloat("CyclicStart"))
        conf_count = Counter_parallel(0)
        s = multiprocessing.Semaphore(1)
        jobs = [multiprocessing.Process(target=worker, name=str(i), args=(s, pool, conf_count, molfile,
                                                                          configfile, testing))
                for i in range(numthreads)]
        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        logger.warning("Normal termination. Generated totally %d conformations." % (conf_count.value()))

    if testing:
        testobj = Testing()
        for molfile in glob.glob(config["Testing"].get("TestMolecules")):
            run_parallel_mc(molfile)
    else:
        run_parallel_mc(molfile)

if __name__ == '__main__':
    print("----------\nUsage: python parallel_mc.py 4 myconfig.cfg\nWhere 4 - the number of threads,\nmyconfig.cfg - "
          "configuration file\n----------")
    if not sys.argv[1].isdigit():
        raise Exception("First argument must represent a number of threads")
    if len(sys.argv) > 2 and sys.argv[2] != "test":
        main(sys.argv[2], int(sys.argv[1]), testing=("test" in sys.argv))
    else:
        main("default.ini", int(sys.argv[1]), testing=("test" in sys.argv))
