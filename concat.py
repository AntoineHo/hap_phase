#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys

from utils import log

def read_mafft(input) :

    haplotypes = {}

    for line in open(input, "r") :
        if line.strip()[0] == ">" :
            name = line.strip().split()[0][1:]
            haplotypes[name] = ""
        else :
            haplotypes[name] += line.strip()

    return haplotypes

def add_missing_samples_and_get_consensus(samples, haplotypes) :

    max_length = max([len(s) for s in haplotypes.values()])
    missing_string = "".join('-' for i in range(max_length))

    seqs = {}

    for sm in samples :
        h1 = sm + "_hap_1"
        h2 = sm + "_hap_2"

        if h1 not in haplotypes.keys() :
            seqs[sm] = missing_string
        else :
            seq = consensus_haplotypes(haplotypes[h1], haplotypes[h2])
            seqs[sm] = seq

    return seqs

def concatenate(samples, all_sequences) :
    seqs = {sm:"" for sm in samples}

    for seqdc in all_sequences :
        for sm, sm_seq in seqdc.items() :
            seqs[sm] += sm_seq

    return seqs

def consensus_haplotypes(h1, h2) :

    consensus = ""
    for i, j in zip(h1, h2) :
        if i == j :
            consensus += i
        elif any(n == "A" for n in [i,j]) and any(n == "T" for n in [i,j]) : # weak
            consensus += "W"
        elif any(n == "C" for n in [i,j]) and any(n == "G" for n in [i,j]) : # strong - CG
            consensus += "S"
        elif any(n == "A" for n in [i,j]) and any(n == "C" for n in [i,j]) : # amino - AC
            consensus += "M"
        elif any(n == "G" for n in [i,j]) and any(n == "T" for n in [i,j]) : # ketone - GT
            consensus += "K"
        elif any(n == "A" for n in [i,j]) and any(n == "G" for n in [i,j]) : # purine - AG
            consensus += "R"
        elif any(n == "C" for n in [i,j]) and any(n == "T" for n in [i,j]) : # pyrimidine - CT
            consensus += "Y"
        else :
            if i == '-' :
                consensus += j # if any haplotype has a gap
            elif j == '-' :
                consensus += i
            else :
                consensus += '-'

    return consensus

def concat(args) :

    # IO
    directory = args.DIR[0]
    sample_dir = args.SAMPLEDIR[0]
    output = args.output[0]

    haplotypes = []

    log("Finding all samples...")
    vcfs = [os.path.join(sample_dir, file) for file in os.listdir(sample_dir) if file.endswith(".phased.vcf")]
    samples = [os.path.basename(fl).split(".")[0] for fl in vcfs]

    log("Finding .mafft files in directory...")
    mafft_files = [os.path.join(directory, fl) for fl in os.listdir(directory) if fl.endswith(".mafft")]

    log("Reading mafft files...")
    all_sequences = []
    for file in mafft_files :
        haplotypes = read_mafft(file)
        seqs = add_missing_samples_and_get_consensus(samples, haplotypes)
        all_sequences.append(seqs)

    log("Concatenating all sequences...")
    final_seqs = concatenate(samples, all_sequences)

    log("Outputting to file...")
    f = open(output, 'w')
    for sm, sequence in final_seqs.items() :
        f.write(">{}\n".format(sm))
        f.write(sequence + "\n")
    f.close()

    log("Done!")
