#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Using phased sequences from same regions in different samples.
Check which regions are phased in how many samples then output fasta for alignments"""

import os
import sys
import argparse

import numpy as np

from utils import log

#### Specific parsers to this step
def parse_bed(bed) :
    regions = []

    for i, line in enumerate(open(bed, "r")) :
        s = line.strip().split("\t")
        name = "{}_{}_{}".format(s[0], int(s[1])+1, s[2]) # vcf coordinates are 1 indexed so need to +1 the bed coordinates
        regions.append(name)
    return regions

def parse_fasta(fasta) :
    """reads a fasta file storing record ids and sequences"""
    seq = {}

    f = open(fasta, "r")
    for line in f :
        if line[0] == ">" :
            name = line[1:].strip().split()[0] # remove the end at first space
            name = "_".join(c for c in name.split("_")[:4])
            seq[name] = ""
        else :
            seq[name] += line.strip()

    f.close()
    return seq


def compare_samples(regions, samples) :
    """For each samples, store record id from fasta and check which are common"""

    #print(regions)
    hap_counts = {k:[] for k in regions} # format = region name : [list of samples phased in this region]
    hap_lengths = {k:[] for k in regions}

    for sample, data in samples.items() :
        dir = data[0]
        h1 = parse_fasta(data[1])
        #h2 = parse_fasta(data[2])

        for name, seq in h1.items() :
            #name = "_".join(c for c in region.split("_")[:4])
            hap_counts[name].append(sample)
            hap_lengths[name].append(len(seq))

    return hap_counts, hap_lengths

def parse_samples(sampledir) :
    """scrape phased fasta files for samples"""

    names = []
    for fl in os.listdir(sampledir) :
        if fl.endswith(".phased.vcf") :
            sample_name = fl.split(".")[0]
            names.append(sample_name)
        else :
            continue

    samples = {}
    for sample in names :
        hap1 = os.path.join(sampledir, sample + ".hap1.fa")
        hap2 = os.path.join(sampledir, sample + ".hap2.fa")
        if not os.path.isfile(hap1) or not os.path.isfile(hap2) :
            raise Exception("Could not find haplotype path for sample {}".format(sample))
        else :
            samples[sample] = (sample_path, hap1, hap2)
    return samples

def write_fasta(filename, samples, reg, sample_list) :

    f = open(filename, "w")
    for sm in sample_list :
        hap1 = parse_fasta(samples[sm][1])
        hap2 = parse_fasta(samples[sm][2])

        f.write(">{}_hap_1\n".format(sm))
        f.write(hap1[reg] + "\n")
        f.write(">{}_hap_2\n".format(sm))
        f.write(hap2[reg] + "\n")

    f.close()
    print("Done region {}".format(reg))

def group(args) :

    fasta = os.path.abspath(args.FASTA[0])
    bed = os.path.abspath(args.BED[0])
    try :
        out = os.path.abspath(args.OUTPUT[0])
        if not os.path.isdir(out) :
            os.makedirs(out)
    except :
        out = os.path.join(os.getcwd(), args.OUTPUT[0])
        out = os.path.abspath(out)
        try :
            os.makedirs(out)
        except :
            raise Exception("Output directory path is invalid")
    #print(args.SAMPLES[0])

    ref = parse_fasta(args.FASTA[0])
    reg = parse_bed(args.BED[0])
    samples = parse_samples(args.SAMPLEDIR[0])
    maxsm = len(samples)
    #print(maxsm)
    #print(reg[0:3])

    # parse samples
    hap_counts, hap_lengths = compare_samples(reg, samples)

    counts = {i:0 for i in range(maxsm)}
    for region in reg :
        sms = hap_counts[region]
        counts[len(sms)] += 1

    print("Histogram: phased regions shared by exactly n samples:")
    for i, c in counts.items() :
        print("{}\t{}".format(i, c))

    print("==================")
    print("Histogram: phased regions shared by at least n samples:")
    for i, c in counts.items() :
        atleast = sum([counts[j] for j in range(i, maxsm)])
        print("{}\t{}".format(i, atleast))

    for reg in hap_counts.keys() :
        concerned_samples = hap_counts[reg]
        concerned_lengths = hap_lengths[reg]

        if len(concerned_samples) >= args.min[0] and np.mean(concerned_lengths) >= args.length :
            filename = reg + ".fasta"
            outfasta = os.path.join(out, filename)
            write_fasta(outfasta, samples, reg, concerned_samples)

    sys.exit(0)
