import configparser, builtins, random, os, traceback, sys, glob
from datetime import date
from shutil import copy2

def writeGeom(config, mol, logger, conf_count):
    logger.info("Conformer was generated successfully.")
    if ".mol" in config["MCDriver"].get("OutputFile"):
        mol.writeToMol(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    elif ".xyz" in config["MCDriver"].get("OutputFile"):
        mol.writeToXyz(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))

def main(configfile, testing = False):
    config = configparser.ConfigParser()
    config.read(configfile)

    builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
    builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
    builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
    from rcrilib import IK_Molecule, SystematicDriver, Testing
    from rcrilib.Helpers import createLogger, log_versions
    do_syssearch_noncyclic = config["SystematicSearch"].getboolean("EnableNonCyclic")
    do_syssearch_cyclic = config["SystematicSearch"].getboolean("EnableCyclic")

    logger = createLogger("MainScript")
    logger.info("Starting Monte-Carlo conformer generation")
    molfile = config["MCDriver"].get("InputFile")
    if not os.path.isfile(molfile):
        raise Exception("Input file \"%s\" is not given" % molfile)

    if testing:
        maxconf = config["Testing"].getint("NumberOfConfs")
    else:
        maxconf = config["MCDriver"].getint("NumberOfConfs")

    def run_recursive_mc(molfile):
        try:
            if testing:
                testobj.reset()
                testobj.check_start()
            log_versions(logger)
            logger.warning("Initializing from file \"%s\"" % molfile)
            mol = IK_Molecule(molfile, config)
            mol.prepare_ik()

            if do_syssearch_cyclic:
                sysdriver_cyc = SystematicDriver(mol.PS.getCyclicDOFs(),
                                                        config["SystematicSearch"].getfloat("CyclicStep"),
                                                        config["SystematicSearch"].getfloat("CyclicStart"))
            if do_syssearch_noncyclic:
                sysdriver_noncyc = SystematicDriver(mol.PS.getNonCyclicDOFs(),
                                                        config["SystematicSearch"].getfloat("NonCyclicStep"),
                                                        config["SystematicSearch"].getfloat("NonCyclicStart"))
            if testing:
                testobj.stop_init()

            conf_count = 0
            while conf_count != maxconf:
                ps = mol.getPS()
                ntry = 0
                logger.warning("Entering refinment stage")
                while ntry < config["MCDriver"].getint("MaxTries") and not ps.success:
                    for item in ps:
                        if item.isContinuous():
                            item.setValue(random.uniform(-3.141592, 3.141592))
                        if item.isDiscrete():
                            item.setValue(None)
                    if ntry == 0 and do_syssearch_cyclic and not sysdriver_cyc.isDone():
                        logger.warning("Making a step in cyclic DOFs")
                        sysdriver_cyc.setValues()
                        sysdriver_cyc.switchToNext()
                    ps = mol.applyPS()
                    ntry += 1
                    if ps.success:
                        if not testing:
                            writeGeom(config, mol, logger, conf_count)
                        else:
                            testobj.stop_good()
                        conf_count += 1
                        break
                    elif testing:
                        testobj.stop_bad()
                    if ps.isEmpty():
                        ntry = config["MCDriver"].getint("MaxTries")
                        break
                if ntry == config["MCDriver"].getint("MaxTries"):
                    mol.perturb_geometry()
                    logger.warning("Couldn't generate the initial conformation. Trying again.")
                    continue
                elif conf_count == maxconf:
                    break

                if do_syssearch_noncyclic:
                    logger.warning("Starting non-cyclic search")
                    while not sysdriver_noncyc.isDone() and conf_count != maxconf:
                        sysdriver_noncyc.setValues()
                        sysdriver_noncyc.switchToNext()
                        ps = mol.applyPS()
                        if ps.success:
                            if not testing:
                                writeGeom(config, mol, logger, conf_count)
                            else:
                                testobj.stop_good()
                            conf_count += 1
                        elif testing:
                            testobj.stop_bad()
                    sysdriver_noncyc.reset()
                if not testing and do_syssearch_cyclic and sysdriver_cyc.isDone() and \
                        config["SystematicSearch"].getboolean("TerminateWhenDone"):
                    logger.info("Systematic search is finished. Terminating.")
                    break
                if not ps.success:
                    logger.info("No luck in this iteration. Trying again.")
                mol.perturb_geometry()
            logger.warning("Normal termination. Generated totally %d conformations." % (conf_count))
            if testing:
                testobj.summarize(molfile)
        except Exception:
            if Exception != KeyboardInterrupt:
                logger.error("Exception. Stack trace:\n" + traceback.format_exc())
                logger.error("Creating bug report.")
                i = 0
                while os.path.isdir("bugreport_%s_%d" % (str(date.today()), i)):
                    i += 1
                bug_dir = "bugreport_%s_%d" % (str(date.today()), i)
                os.mkdir(bug_dir)
                try:
                    copy2(molfile, bug_dir + "/starting_geom.mol")
                    copy2(config["MCDriver"].get("LogFile"), bug_dir + "/logfile")
                except:
                    pass

                # Write values of all the DOFs
                torfile = open(bug_dir + "/torvalues.csv", "w")
                ps = mol.getPS()
                torlines = ["atom1,atom2,side1,side2,value"]
                for param in ps:
                    try:
                        torlines.append("%d,%d,%d,%d,%f" % (param.atoms[0], param.atoms[1],
                                                            param.sides[0], param.sides[1],
                                                            param.value))
                    except:
                        torlines.append("Error accessing parameter attributes")
                torfile.write("\n".join(torlines))
                torfile.close()
        finally:
            if config["MCDriver"].getboolean("AutoCleanup") and os.path.isfile(config["MCDriver"].get("LogFile")):
                os.remove(config["MCDriver"].get("LogFile"))

    if testing:
        testobj = Testing()
        for molfile in glob.glob(config["Testing"].get("TestMolecules")):
            run_recursive_mc(molfile)
    else:
        run_recursive_mc(molfile)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] != "test":
        main(sys.argv[1], testing=("test" in sys.argv))
    else:
        main("default.ini", testing=("test" in sys.argv))
