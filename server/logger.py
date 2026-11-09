# -*- coding: utf-8 -*-
"""
Shared logging setup for the hvFinder server scripts.
"""
import logging
import sys

def setup_logger(name, log_file, level=logging.INFO):
    """
    Configure the common log format and level.
    """
    formatter = logging.Formatter('<%(asctime)s> [%(name)s] %(levelname)s: %(message)s', 
                                  datefmt='%Y-%m-%d %H:%M:%S')

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.addHandler(file_handler)
    
    return logger