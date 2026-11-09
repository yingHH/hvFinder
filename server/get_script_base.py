#!/usr/bin/env python3
"""
Print the directory this script lives in, so clients can locate the
installation.
"""

import os

def get_script_base():
    """
    Return the directory of this script.
    """
    return os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    print(get_script_base())
