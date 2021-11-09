import random, os, multiprocessing, configparser, builtins
from copy import deepcopy

from rcrilib import IK_Molecule, SystematicDriver, SystematicDriver_parallel, Confcounter_parallel
from rcrilib.Helpers import log_versions, cause, Confcounter

def writeGeom(config, mol, logger, conf_count):
    logger.info("Conformer was generated successfully.")
    if ".mol" in config["MCDriver"].get("OutputFile"):
        mol.writeToMol(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    elif ".xyz" in config["MCDriver"].get("OutputFile"):
        mol.writeToXyz(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    # logger.info("Excel: %s created" % config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    # logger.info("Short: %s created" % config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))

def run_parallel_mc(molfile, config, configfile, numthreads, logger, timer, testing):
    if not os.path.isfile(molfile):
        raise Exception("Input file \"%s\" is not found" % molfile)

    logger.error("Starting parallel CS for \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    pool = SystematicDriver_parallel(mol.PS.getCyclicDOFs(),
                                     config["SystematicSearch"].getfloat("CyclicStep"),
                                     config["SystematicSearch"].getfloat("CyclicStart"))
    if testing:
        maxconf = config["Testing"].getint("NumberOfConfs")
    else:
        maxconf = config["MCDriver"].getint("NumberOfConfs")
    conf_count = Confcounter_parallel(maxconf)

    with multiprocessing.Manager() as manager:
        timer.reset(manager)
        timer.check_start()
        jobs = [multiprocessing.Process(target=parallel_mc_worker, name=str(i), args=(pool, conf_count, molfile,
                                                                          configfile, timer, testing))
                for i in range(numthreads)]
        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        timer.format_lists()
    logger.warning("Normal termination. Generated totally %d conformations." % (conf_count.value()))


def parallel_mc_worker(pool, conf_count, molfile, configfile, timer, testing):
    config = configparser.ConfigParser()
    config.read(configfile)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    builtins.scounter_linecount = 0

    from rcrilib import IK_Molecule, SystematicDriver
    from rcrilib.Helpers import createLogger, log_versions
    do_syssearch_noncyclic = config["SystematicSearch"].getboolean("EnableNonCyclic")
    do_syssearch_cyclic = config["SystematicSearch"].getboolean("EnableCyclic")

    logger = createLogger("WorkerScript")
    logger.warning("Starting IK-MC conformer generation with bruteforce of DDOFs")
    if not os.path.isfile(molfile):
        raise Exception("Input file \"%s\" is not found" % molfile)

    log_versions(logger)
    logger.warning("Initializing from file \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    if do_syssearch_cyclic:
        sysdriver_cyc = pool
    else:
        sysdriver_cyc = None

    if do_syssearch_noncyclic:
        sysdriver_noncyc = SystematicDriver(mol.PS.getNonCyclicDOFs(),
                                            config["SystematicSearch"].getfloat("NonCyclicStep"),
                                            config["SystematicSearch"].getfloat("NonCyclicStart"))
    else:
        sysdriver_noncyc = None

    timer.stop_init()
    run_ddof_mc(mol, config, logger, timer, conf_count, testing, sysdriver_cyc, sysdriver_noncyc, True)

def singthr_mc(molfile, config, logger, timer, testing):
    timer.check_start()
    do_syssearch_noncyclic = config["SystematicSearch"].getboolean("EnableNonCyclic")
    do_syssearch_cyclic = config["SystematicSearch"].getboolean("EnableCyclic")

    if testing:
        maxconf = config["Testing"].getint("NumberOfConfs")
    else:
        maxconf = config["MCDriver"].getint("NumberOfConfs")
    conf_count = Confcounter(maxconf)

    log_versions(logger)
    logger.error("Initializing from file \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    if do_syssearch_cyclic:
        sysdriver_cyc = SystematicDriver(mol.PS.getCyclicDOFs(),
                                         config["SystematicSearch"].getfloat("CyclicStep"),
                                         config["SystematicSearch"].getfloat("CyclicStart"))
    else:
        sysdriver_cyc = None

    if do_syssearch_noncyclic:
        sysdriver_noncyc = SystematicDriver(mol.PS.getNonCyclicDOFs(),
                                            config["SystematicSearch"].getfloat("NonCyclicStep"),
                                            config["SystematicSearch"].getfloat("NonCyclicStart"))
    else:
        sysdriver_noncyc = None

    timer.stop_init()
    run_ddof_mc(mol, config, logger, timer, conf_count, testing, sysdriver_cyc, sysdriver_noncyc, False)

def run_ddof_mc(mol, config, logger, timer, conf_count, testing, sysdriver_cyc, sysdriver_noncyc, parallel):
    do_syssearch_cyclic = sysdriver_cyc is not None
    do_syssearch_noncyclic = sysdriver_noncyc is not None

    major_iterations = 0
    if parallel and do_syssearch_cyclic:
        cyclicdofs = mol.PS.getCyclicDOFs()

    while not conf_count.isFinished() and major_iterations != config["MCDriver"].getint("MaxIter"):
        logger.error("Generated %d conformations. Goal = %d" % (conf_count, conf_count.maxconf))
        major_iterations += 1
        ps = mol.getPS()
        ps.success = False
        mol.discr_cp = -1
        it = 0

        # with open("randomstate.pickle", "wb") as f:
        #     pickle.dump(random.getstate(), f)
        # with open("randomstate.pickle", "rb") as f:
        #     random.setstate(pickle.load(f))

        logger.warning("Generating the initial conformation")

        while it != config["MCDriver"].getint("MaxTries") and not ps.success:
            for item in ps:
                if item.isContinuous():
                    item.setValue(random.uniform(-3.141592, 3.141592))
                if item.isDiscrete():
                    item.setValue(None)
            if it == 0 and do_syssearch_cyclic and not sysdriver_cyc.isDone():
                logger.warning("Making a step in cyclic DOFs")
                if parallel:
                    sysdriver_cyc.setValues(cyclicdofs)
                else:
                    sysdriver_cyc.setValues()
                sysdriver_cyc.switchToNext()
            ps = mol.applyPS()
            it += 1
            if ps.success:
                timer.stop_good()
                conf_count.increment()
                if not testing:
                    writeGeom(config, mol, logger, conf_count)
                break
            else:
                timer.stop_bad()
                if ps.isEmpty():
                    it = config["MCDriver"].getint("MaxTries")
                    break
        if it == config["MCDriver"].getint("MaxTries"):
            mol.perturb_geometry()
            logger.warning("Couldn't generate the initial conformation. Trying again.")
            continue
        elif conf_count.isFinished():
            break

        backup_ps = deepcopy(mol.getPS())
        if do_syssearch_noncyclic:
            logger.warning("Starting non-cyclic search")
            while not sysdriver_noncyc.isDone() and not conf_count.isFinished():
                sysdriver_noncyc.setValues()
                sysdriver_noncyc.switchToNext()
                ps = mol.applyPS()
                if ps.success:
                    timer.stop_good()
                    conf_count.increment()
                    if not testing:
                        writeGeom(config, mol, logger, conf_count)
                else:
                    timer.stop_bad()
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
        while not mol.done_all_discr() and not conf_count.isFinished():
            # logger.warning("New DDOF combination")
            # logger.info("Excel: NEW APS(ID)")
            # logger.info("Short: NEW APS(ID)")
            # logger.info("Excel: %s" % mol.ccounter.logstate())
            # check = "Short: %s" % mol.ccounter.logstate(nocount=True)
            # check2 = mol.ccounter.logstate()
            # if "line #2144" in check2:
            #     logger.info("Here I am")
            # logger.info("Short: %s" % check2)
            check = mol.ccounter.logstate(nocount=True)
            ps = mol.applyPS(increase_discrete=True)
            # logger.info("Excel: %s" % mol.ccounter.logstate())
            # logger.info("Short: %s" % mol.ccounter.logstate())
            mol.ccounter.backup_ddof()
            # logger.info("Excel: --- %s;%s\n" % (repr(ps.cause), repr(ps.success)))
            # logger.info("Short: --- %s;%s\n" % (repr(ps.cause), repr(ps.success)))
            mol.ccounter.checkState()
            if prev == check:
                raise Exception("DDOF state repetition!!!")
            prev = check
            trycount += 1
            it = 0

            if ps.success:
                timer.stop_good()
                conf_count.increment()
                if not testing:
                    writeGeom(config, mol, logger, conf_count)

            while it != config["MCDriver"].getint("MaxTries") and not ps.success and \
                    ps.cause != cause.geomoverlap_ddof and ps.cause != cause.zerosolutions_ddof:
                for item in ps:
                    if item.isContinuous():
                        logger.info("Setting random torsion on " + repr(item.atoms))
                        item.setValue(random.uniform(-3.141592, 3.141592))
                # logger.info("Excel: New try:")
                # check3 = mol.ccounter.logstate()
                # if "line #3348" in check3:
                #     logger.info("Here I am")
                # logger.info("Excel: %s" % check3)
                ps = mol.applyPS()
                # logger.info("Excel: %s" % mol.ccounter.logstate())
                # logger.info("Excel: --- %s;%s\n" % (repr(ps.cause), repr(ps.success)))
                it += 1
                if not ps.success and ps.isEmpty():
                    break
                if ps.success:
                    timer.stop_good()
                    conf_count.increment()
                    if not testing:
                        writeGeom(config, mol, logger, conf_count)
                    break
                else:
                    timer.stop_bad()

            if do_syssearch_noncyclic:
                logger.warning("Starting non-cyclic search")
                while not sysdriver_noncyc.isDone() and not conf_count.isFinished():
                    sysdriver_noncyc.setValues()
                    sysdriver_noncyc.switchToNext()
                    ps = mol.applyPS()
                    if ps.success:
                        timer.stop_good()
                        conf_count.increment()
                        if not testing:
                            writeGeom(config, mol, logger, conf_count)
                    else:
                        timer.stop_bad()
                sysdriver_noncyc.reset()

            ps = mol.getPS()
            for i in range(len(ps)):
                if ps[i].isContinuous():
                    if ps[i].value != backup_ps[i].value:
                        ps[i].setValue(backup_ps[i].value)
            mol.ccounter.cleanup_temporary_ddofs()
            mol.ccounter.recordState()
        logger.warning("DDOF bruteforce is finished")
        mol.perturb_geometry()
        if not testing and do_syssearch_cyclic and sysdriver_cyc.isDone() and \
                config["SystematicSearch"].getboolean("TerminateWhenDone"):
            logger.warning("Systematic search is finished. Terminating.")
            break
    logger.warning("Normal termination. Generated totally %d conformations." % (conf_count))
