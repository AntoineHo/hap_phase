#!/usr/bin/python
# -*- coding: utf-8 -*-

# Standard library
import os
import sys
import subprocess

# Modules
from utils import log, run, which
from parsers import parse_samples

##################### MAIN #####################

def run_sambamba(bam, outdir, threads) :

    bam_basename = os.path.basename(bam)
    output = os.path.join(outdir, bam_basename + ".cov")
    comp = output + ".gz"

    if not os.path.exists(comp) :
        # Run sambamba if not run already
        if not os.path.exists(output) :
            cmd = "sambamba depth base -t {threads} --min-coverage=0 --min-base-quality=0 {bam} > {output}" #" | gzip --best > {output}"
            dc_args = {"threads":threads, "bam":bam, "output":output}
            cmd = cmd.format(**dc_args)
            print(cmd)
            run(cmd)

        # Run gzip after sambamba
        cmd = "gzip --best {input} > {output}" #" | gzip --best > {output}"
        dc_args = {"input":output, "output":comp}
        cmd = cmd.format(**dc_args)
        run(cmd)
        
    else :
        log("Found coverage file at: {}".format(comp), ret = False)

    return comp

def get_coverage(args) :
    """Runs hapcut2: require bcftools, HAPCUT2 and extractHAIRS in $PATH"""

    # IO
    sample_file = args.SAMPLES[0]
    outdir = args.outdir[0]
    threads = args.threads[0]

    if not os.path.exists(outdir) :
        os.makedirs(outdir) # Create directory following path
    else :
        log("WARNING: Output directory already exists: {}.".format(outdir))

    print("# Sambamba runner tool")
    print("Input sample list:\t{0}".format(sample_file))
    print("Output path:\t\t{0}".format(outdir))
    print("Other parameters:")
    print("\t- Threads: {}".format(threads))
    print(
        "===============================================================================\n"
    )

    log("Parsing samples file")
    samples = parse_samples(sample_file)

    #log("Parsing .bed file")
    #regions = parse_bed(bed_file)

    log("Running sambamba depth on all samples...")
    i = 1
    for sample, bam in samples.items() :
        log("Sample {}/{}: {}".format(i, len(samples), sample))
        run_sambamba(bam, outdir, threads)
        i += 1

    log("Done!")
