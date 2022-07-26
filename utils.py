#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess
from time import localtime, strftime

def log(string: str):
    """Prints a log"""
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))

def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()
    return proc.returncode

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return False
