# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2026-01-21 15:55:16
Last Modified by: Ying Huang
Last Modified time: 2026-01-21 15:55:16

Description: 
统一日志记录模块
"""
import logging
import sys
import os

def setup_logger(name, log_file, level=logging.INFO):
    """设置统一的日志格式"""
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