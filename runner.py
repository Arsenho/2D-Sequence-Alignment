#!/usr/bin/python3

import re
import imp
import os


def filter_shfile(name):
    return re.match("^\w+.sh$", name)


hdir = ["matrics/sequential", "matrics/parallel"]

Hierarchies = {}
cmd = "cp -r {}/* ."

for d in hdir:
    os.system(cmd.format(d))

for d in hdir:
    files = [f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]
    # py_files = filter(filter_pyfile, files)
    sh_files = filter(filter_shfile, files)
    # py_files += [file for file in sh_files]

    for f in sh_files:
        cmd = "sbatch {}"
        os.system(cmd.format(f))
