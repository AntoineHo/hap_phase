#!/usr/bin/python
# -*- coding: utf-8 -*-

# Standard library
import os
import sys
import subprocess

# Modules
from utils import log, run, which
#from parsers import parse_bed

##################### PARSING FILES #####################

def parse_samples(sample_file) :
    """Parse a sample file like: sample_name\tbam_file_path"""
    samples = {}

    f = open(sample_file, "r")
    for line in f :
        if line.startswith("#") :
            continue

        s = line.strip().split("\t")
        sm = s[0]
        bam = s[1]

        if not os.path.exists(bam) :
            log("ERROR: Could not find .bam filepath: {}".format(bam))
            sys.exit(1)

        samples[sm] = bam

    f.close()

    return samples

##################### MAIN #####################

def filter_input_vcf(sample, bed, vcf, output, threads) :

    output = os.path.join(outdir, sample + ".for_phasing.vcf")
    print(output)
    cmd = "bcftools view -Ov -R {bed} -s {sample} --threads {threads} -g het --exclude-uncalled -m2 -M2 -o {output} {vcf}"
    dc_args = {"bed":bed, "vcf":vcf, "sample":sample, "threads":threads, "output":output}
    cmd = cmd.format(**dc_args)
    print(cmd)
    #run(cmd)

    return output

def run_extractHAIRS(bam) :
    pass

def run_HAPCUT2() :
    pass

def phase_vcf(args) :
    """Runs hapcut2: require bcftools, HAPCUT2 and extractHAIRS in $PATH"""

    # IO
    sample_file = args.SAMPLES[0]
    vcf_file = args.VCF[0]
    bed_file = args.BED[0]
    outdir = args.outdir[0]
    threads = args.threads[0]

    print("# Hapcut2 runner tool")
    print("Input VCF to phase:\t{0}".format(vcf))
    print("Input BED with regions to phase:\t{0}".format(vcf))
    print("Input sample list:\t{0}".format(samples))
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

    log("Running bcftools, extractHAIRS and HAPCUT2 on all samples...")
    for sm, bam in samples.items() :
        log("Starting phasing sample: {}".format(sm))


        filtered_vcf_path = filter_input_vcf(sample, bed_file, vcf_file, output, threads)

        log("Finished phasing sample: {}".format(sm), ret = False) # No return char (line sticks to previous log/print line)



    log("Done!")
