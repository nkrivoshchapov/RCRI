import configparser
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

from rcrilib import IK_Molecule
from rcrilib.Helpers import createLogger, log_versions

logger = createLogger("MainScript")
logger.info("Starting Monte-Carlo conformer generation")
deg2rad = 0.0174532925199432957692
rad2deg = 57.295779513082320877

molfile = config["MCDriver"].get("InputFile")
if not os.path.isfile(molfile):
    raise Exception("Input file \"%s\" is not given" % molfile)

def writeGeom():
    logger.info("Conformer was generated successfully.")
    if ".mol" in config["MCDriver"].get("OutputFile"):
        mol.writeToMol(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))
    elif ".xyz" in config["MCDriver"].get("OutputFile"):
        mol.writeToXyz(config["MCDriver"].get("OutputFile").replace("%d", str(conf_count)))

try:
    log_versions(logger)
    logger.info("Initializing from file \"%s\"" % molfile)
    mol = IK_Molecule(molfile, config)
    mol.prepare_ik()
    conf_count = 0
    while conf_count < config["MCDriver"].getint("NumberOfConfs"):
        niter = 0
        while niter != config["MCDriver"].getint("MaxIter"):
            ps = mol.getPS()
            ntry = 0
            while ntry < config["MCDriver"].getint("MaxTries") and not ps.isEmpty():
                for item in ps:
                    if item.isContinuous():
                        item.setValue(random.uniform(-3.141592, 3.141592))
                    if item.isDiscrete():
                        item.setValue(None)
                ps = mol.applyPS()
                ntry += 1
            if ps.success:
                writeGeom()
                conf_count += 1
                break
            logger.info("No luck in this iteration. Trying again.")
            mol.perturb_geometry()
            niter += 1
        if niter == config["MCDriver"].getint("MaxIter"):
            logger.warning("Reached MaxIter")
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
