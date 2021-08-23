import configparser,pickle
from copy import deepcopy
from datetime import date
import random, os, traceback
from shutil import copy2

config = configparser.ConfigParser()
if len(sys.argv) > 1:
    config.read(sys.argv[1])
else:
    config.read("default.ini")

import builtins
builtins.rcrilib_pl_console = config["MCDriver"].get("PrintLevelConsole")
builtins.rcrilib_pl_file = config["MCDriver"].get("PrintLevelFile")
builtins.rcrilib_logfilename = config["MCDriver"].get("LogFile")
builtins.scounter_linecount = 0

from rcrilib import IK_Molecule
from rcrilib.Helpers import createLogger, log_versions, cause

logger = createLogger("MainScript")
logger.warning("Starting IK-MC conformer generation with bruteforce of DDOFs")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

molfile = config["MCDriver"].get("InputFile")
if not os.path.isfile(molfile):
    raise Exception("Input file \"%s\" is not found" % molfile)

def writeGeom():
    if ".mol" in config["MCDriver"].get("OutputFile"):
        mol.writeToMol(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
        logger.warning("Writing file " + config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    elif ".xyz" in config["MCDriver"].get("OutputFile"):
        mol.writeToXyz(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
        logger.warning("Writing file " + config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    logger.info("Conformer was generated successfully.")

try:
    with open("randomstate.pickle", "wb") as f:
        pickle.dump(random.getstate(), f)
    log_versions(logger)
    logger.info("Initializing from file \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    conf_count = 0
    while conf_count < config["MCDriver"].getint("NumberOfConfs"):
        logger.warning("Generated %d conformations. Goal = %d" % (conf_count,
                                                                  config["MCDriver"].getint("NumberOfConfs")))
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
            ps = mol.applyPS()
            it += 1
            if ps.success:
                writeGeom()
                conf_count += 1
                break
            elif ps.isEmpty():
                it = config["MCDriver"].getint("MaxTries")
                break
        if it == config["MCDriver"].getint("MaxTries"):
            mol.perturb_geometry()
            logger.warning("Couldn't generate the initial conformation. Trying again.")
            continue
        backup_ps = deepcopy(mol.getPS())
        mol.startDiscreteRun()
        trycount = 0
        prev = ""
        logger.warning("Starting bruteforce of DDOF combinations")
        while not mol.done_all_discr():
            logger.warning("New DDOF combination")
            check = mol.ccounter.logstate(nocount=True)
            ps = mol.applyPS(increase_discrete = True)
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
                        logger.info("Setting random torsion on "+repr(item.atoms))
                        item.setValue(random.uniform(-3.141592, 3.141592))
                ps = mol.applyPS()
                it += 1
                if not ps.success and ps.isEmpty():
                    break
            if ps.success:
                writeGeom()
                conf_count += 1
            ps = mol.getPS()
            for i in range(len(ps)):
                if ps[i].isContinuous():
                    if ps[i].value != backup_ps[i].value:
                        ps[i].setValue(backup_ps[i].value)
            mol.ccounter.recordState()
        logger.warning("DDOF bruteforce is finished")
        mol.perturb_geometry()
    logger.warning("Normal termination. Generated totally %d conformations." % (conf_count))
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
            copy2("randomstate.pickle", bug_dir + "/randomstate.pickle")
            copy2(__file__, bug_dir + "/driver.py")
            copy2(molfile, bug_dir + "/starting_geom.mol")
            copy2(config["MCDriver"].get("LogFile"), bug_dir + "/logfile")
        except:
            pass

        # Write values of all the DOFs
        torfile = open(bug_dir + "/parameters_continuous.csv", "w")
        ps = mol.getPS()
        torlines = ["atom1,atom2,side1,side2,value"]
        for param in ps:
            try:
                if param.isContinuous():
                    torlines.append("%d,%d,%d,%d,%f" % (param.atoms[0], param.atoms[1],
                                                        param.sides[0], param.sides[1],
                                                        param.value))
            except:
                torlines.append("Error accessing parameter attributes")
        torfile.close()
finally:
    if config["MCDriver"].getboolean("AutoCleanup") and os.path.isfile(config["MCDriver"].get("LogFile")):
        os.remove(config["MCDriver"].get("LogFile"))
