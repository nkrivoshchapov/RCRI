import configparser, builtins, random, os, traceback, sys, time, glob
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
    from rcrilib import IK_Molecule, SystematicDriver, JointDriver, Testing
    from rcrilib.Helpers import createLogger, log_versions
    do_syssearch_noncyclic = config["SystematicSearch"].getboolean("EnableNonCyclic")
    do_syssearch_cyclic = config["SystematicSearch"].getboolean("EnableCyclic")

    logger = createLogger("MainScript")
    logger.info("Starting Monte-Carlo conformer generation")
    if not testing:
        molfile = config["MCDriver"].get("InputFile")
        if not os.path.isfile(molfile):
            raise Exception("Input file \"%s\" is not given" % molfile)

    if testing:
        maxconf = config["Testing"].getint("NumberOfConfs")
    else:
        maxconf = config["MCDriver"].getint("NumberOfConfs")

    def run_simple_mc(molfile):
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
            if do_syssearch_noncyclic and do_syssearch_cyclic:
                sysdriver = JointDriver(sysdriver_cyc, sysdriver_noncyc)
            elif do_syssearch_noncyclic:
                sysdriver = sysdriver_noncyc
            elif do_syssearch_cyclic:
                sysdriver = sysdriver_cyc
            else:
                sysdriver = None
            if testing:
                testobj.stop_init()

            conf_count = 0
            while conf_count != maxconf:
                ps = mol.getPS()
                for item in ps:
                    if item.isContinuous():
                        item.setValue(random.uniform(-3.141592, 3.141592))
                    if item.isDiscrete():
                        item.setValue(None)
                if sysdriver is not None and not sysdriver.isDone():
                    sysdriver.setValues()
                    sysdriver.switchToNext()
                ps = mol.applyPS()
                if ps.success:
                    if not testing:
                        writeGeom(config, mol, logger, conf_count)
                    else:
                        testobj.stop_good()
                    conf_count += 1
                    mol.perturb_geometry()
                    continue
                elif testing:
                    testobj.stop_bad()

                if not testing and sysdriver is not None and sysdriver.isDone() and \
                                                            config["SystematicSearch"].getboolean("TerminateWhenDone"):
                    logger.info("Systematic search is finished. Terminating.")
                    break
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
            run_simple_mc(molfile)
    else:
        run_simple_mc(molfile)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] != "test":
        main(sys.argv[1], testing=("test" in sys.argv))
    else:
        main("default.ini", testing=("test" in sys.argv))
