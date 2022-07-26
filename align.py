#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Using grouped by region fasta file, run mafft"""

import os
import sys

from utils import run

def align(args) :

    if not os.path.isdir(args.OUTPUT[0]) :
        os.makedirs(args.OUTPUT[0])

    outputdir = os.path.abspath(args.OUTPUT[0])

    if os.path.isdir(os.path.abspath(args.DIR[0])) :
        inputdir = os.path.abspath(args.DIR[0])
    else :
        sys.exit(1)

    for region_file in os.listdir(inputdir) :
        input_file_path = os.path.join(inputdir, region_file)
        if os.path.isfile(input_file_path) and input_file_path.endswith(".fasta") :
            bn = os.path.splitext(os.path.basename(input_file_path))[0]
            output_file_path = os.path.join(outputdir, bn+".mafft")
            # run mafft
            cmd = "mafft --thread 4 {} > {}".format(input_file_path, output_file_path)
            run(cmd)

    sys.exit(0)
