#!/usr/bin/env python3
"""
返回当前脚本所在的目录路径
用于动态获取脚本目录，以查找默认数据文件
"""

import os

def get_script_base():
    """返回脚本所在目录的路径"""
    return os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    print(get_script_base())
