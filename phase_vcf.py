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

def filter_input_vcf(sample, bed, vcf, threads, outdir) :

    output = os.path.join(outdir, sample + ".for_phasing.vcf")
    if os.path.exists(output) :
        return output

    print(output)
    cmd = "bcftools view -Ov -R {bed} -s {sample} --threads {threads} -g het --exclude-uncalled -m2 -M2 -o {output} {vcf}"
    dc_args = {"bed":bed, "vcf":vcf, "sample":sample, "threads":threads, "output":output}
    cmd = cmd.format(**dc_args)
    print(cmd)
    #run(cmd)

    return output

def run_extractHAIRS(bam, filtered_vcf, outdir) :

    output = os.path.join(outdir, sample + ".fragments")
    if os.path.exists(output) :
        return output

    cmd = "extractHAIRS --bam {bam} --VCF {vcf} --out {fragments}"
    dc_args = {"vcf":filtered_vcf, "bam":bam, "fragments":output}
    cmd = cmd.format(**dc_args)
    print(cmd)
    #run(cmd)

    return output

def run_HAPCUT2(fragments, filtered_vcf, outdir) :

    output = os.path.join(outdir, sample + ".hapcut2")
    if os.path.exists(output) :
        return output

    cmd = "HAPCUT2 --fragments {fragments} --VCF {vcf} --out {output} --outvcf 1"
    dc_args = {"vcf":filtered_vcf, "fragments":fragments, "output":output}
    cmd = cmd.format(**dc_args)
    print(cmd)
    #run(cmd)

def phase_vcf(args) :
    """Runs hapcut2: require bcftools, HAPCUT2 and extractHAIRS in $PATH"""

    # IO
    sample_file = args.SAMPLES[0]
    vcf_file = args.VCF[0]
    bed_file = args.BED[0]
    outdir = args.outdir[0]
    threads = args.threads[0]

    if not os.path.exists(outdir) :
        os.makedirs(outdir) # Create directory following path
    else :
        log("WARNING: Output directory already exists: {}.".format(outdir))

    print("# Hapcut2 runner tool")
    print("Input VCF to phase:\t{0}".format(vcf_file))
    print("Input BED with regions to phase:\t{0}".format(bed_file))
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

    log("Running bcftools, extractHAIRS and HAPCUT2 on all samples...")
    i = 1
    for sample, bam in samples.items() :

        log("Sample {}/{}: {}".format(i, len(samples), sample))
        log("Phasing vcf...", ret = False) # No return char (line sticks to previous log/print line)
        filtered_vcf = filter_input_vcf(sample, bed_file, vcf_file, threads, outdir)
        log("Finished phasing vcf", ret = False)
        log("Running extractHAIRS", ret = False)
        fragments = run_extractHAIRS(bam, filtered_vcf, outdir)
        log("Running HAPCUT2", ret = False)
        run_HAPCUT2(fragments, filtered_vcf, outdir)
        log("Finished HAPCUT2", ret = False)

        i += 1

    log("Done!")
