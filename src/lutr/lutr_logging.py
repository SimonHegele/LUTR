"""
Module Name:    setup_logging.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3

Description:    Provides functions:
                - logging_setup()   configures logging
                - logged_worker()   wrapper function to log exceptions from subprocesses 
"""

from logging         import CRITICAL, DEBUG, ERROR, INFO, WARNING, FileHandler, StreamHandler, basicConfig, critical
from os              import path
from sys             import stdout

def logging_setup(args):

    match args.log_level:
        case "debug":
            level=DEBUG
        case "info":
            level=INFO
        case "warning":
            level=WARNING
        case "error":
            level=ERROR
        case "critical":
            level=CRITICAL
    
    file_handler   = FileHandler(filename=path.join(args.outdir, "lutr.log"))
    stdout_handler = StreamHandler(stream=stdout)

    basicConfig(level    = level,
                format   = "%(asctime)s %(levelname)s [%(module)s] %(message)s",
                datefmt  = "%d-%m-%Y %H:%M:%S",
                handlers=[file_handler, stdout_handler])
    
def log_wrapper(func, args):
    
    try:
        if len(args) > 1:
            return func(*args)
        else:
            return func(args)
    except Exception as e:
        critical(e)
        raise e