import logging, numpy, networkx, scipy, builtins

def str_to_pl(text):
    if text.lower() == "debug":
        pl = logging.DEBUG
    elif text.lower() == "info":
        pl = logging.INFO
    elif text.lower() == "warning":
        pl = logging.WARNING
    elif text.lower() == "error":
        pl = logging.ERROR
    elif text.lower() == "none":
        pl = 0
    else:
        raise Exception("Print level \"%s\" is unknown" % text)
    return pl

def createLogger(classname):
    pl_console = logging.INFO
    pl_file = logging.INFO
    min_pl = logging.INFO

    if hasattr(builtins, "rcrilib_pl_console"):
        plc_text = builtins.rcrilib_pl_console.replace(" ", "").lower()
        pl_console = str_to_pl(plc_text)
        if min_pl > pl_console and pl_console != 0:
            min_pl = pl_console

    if hasattr(builtins, "rcrilib_pl_file"):
        plf_text = builtins.rcrilib_pl_file.replace(" ", "").lower()
        pl_file = str_to_pl(plf_text)
        if min_pl > pl_file and pl_file != 0:
            min_pl = pl_file

    logger = logging.getLogger(classname)
    logger.setLevel(min_pl)
    formatter = logging.Formatter("%(name)s:%(levelname)s  %(message)s")

    if pl_file != 0:
        try:
            file_handler = logging.FileHandler(builtins.rcrilib_logfilename)
        except:
            file_handler = logging.FileHandler("rcrilib_run.log")
        file_handler.setLevel(pl_file)
        file_handler.setFormatter(formatter)

    if pl_console != 0:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(pl_console)
        stream_handler.setFormatter(formatter)

    if pl_file != 0:
        logger.addHandler(file_handler)
    if pl_console != 0:
        logger.addHandler(stream_handler)
    return logger

def log_versions(logger):
    logger.debug("-" * 20)
    logger.debug("Package vesions:")
    logger.debug("Networkx: " + repr(networkx.__version__))
    logger.debug("Numpy: " + repr(numpy.__version__))
    logger.debug("Scipy: " + repr(scipy.__version__))
    logger.debug("-" * 20)
