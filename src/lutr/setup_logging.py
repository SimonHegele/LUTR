"""
Module Name:    setup_logging.py
Author:         Simon Hegele
Date:           2025-07-24
Version:        1.0
License:        GPL-3

Description:    Provides functions:
                - logging_setup() to configure logging
"""

import logging
import os
import sys

def logging_setup(args):

    match args.log_level:
        case "debug":
            level=logging.DEBUG
        case "info":
            level=logging.INFO
        case "warning":
            level=logging.WARNING
        case "error":
            level=logging.ERROR
        case "critical":
            level=logging.CRITICAL
    
    file_handler   = logging.FileHandler(filename=os.path.join(args.outdir, "lutr.log"))
    stdout_handler = logging.StreamHandler(stream=sys.stdout)

    logging.basicConfig(level    = level,
                        format   = "%(asctime)s %(levelname)s [%(module)s] %(message)s",
                        datefmt  = "%d-%m-%Y %H:%M:%S",
                        handlers=[file_handler, stdout_handler]
                        )
