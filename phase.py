#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Using HapCut2 phased blocks and VCF and an haploid assembly. Phases regions from a .bed file. Filtering regions based on criteria such as coverage, AF, mismatch quality, etc. Requires an all-site .cov file from the haploid assembly"""

#import errno
import os
import sys

import pandas as pd
import numpy as np
from pysam import VariantFile

from utils import log, run, which
from parsers import parse_bed, parse_fasta_generic

##################### PARSING FILES #####################

def parse_samples(phasedir) :
    # format: {sample: {vcf:vcf, blocks:blocks}}

    files = [os.path.join(phasedir, file) for file in os.listdir(phasedir)]

    # Fetch sample names
    sample_names = []
    for fl in files :
        if fl.endswith(".fragments") :
            sample = fl.split(".")[0]
            sample_names.append(sample)
        else :
            continue

    # Fetch sample files
    samples = {sm:{} for sm in sample_names}
    for fl in files :
        if fl.endswith(".hapcut2") :
            sample = fl.split(".")[0]
            samples[sample]["blocks"] = fl
        elif fl.endswith(".hapcut2.phased.VCF") :
            sample = fl.split(".")[0]
            samples[sample]["vcf"] = fl
        else :
            continue

    return samples

def parse_haplotype_blocks(phased_blocks) :
    """Reads the haplotype blocks informations from HapCut2 output"""
    data = {"bid":[],"vid":[],"al1":[],"al2":[],"chr":[],"pos":[],"ref":[],"alt":[],"gen":[],"mis":[],"inforeads":[]}

    bid = 0
    for n, line in enumerate(open(phased_blocks, "r")) :
        if line[0] == "*" :
            continue
        elif line[:5] == "BLOCK" :
            bid += 1
        else :
            s = line.strip().split()

            data["bid"].append(bid)
            data["vid"].append(int(s[0]))
            if s[1] != "-"  :
                data["al1"].append(int(s[1]))
            else :
                data["al1"].append(-1)

            if s[2] != "-" :
                data["al2"].append(int(s[2]))
            else :
                data["al2"].append(-1)

            data["chr"].append(s[3])
            data["pos"].append(int(s[4]))
            data["ref"].append(s[5])
            data["alt"].append(s[6])
            data["gen"].append(s[7])
            data["mis"].append(float(s[10]))
            data["inforeads"].append(int(s[11]))

    df = pd.DataFrame().from_dict(data)

    haplotypes = []
    for n, hap in df.groupby("bid").agg({"chr":"first", "pos":["min", "max"], "bid":"first"}).iterrows() :
        haplotypes.append( (hap["chr", "first"], hap["pos", "min"], hap["pos", "max"], hap["bid", "first"]) )

    return df, haplotypes

def get_overlap(x1, x2, y1, y2, ) :
    """From two sets of 1d coordinates, returns the amount of overlap"""
    if y2 >= x1 and y1 <= x2 :
        return max(0, min(x2,y2) - max(x1,y1))
    else :
        return 0

def find_overlaps(regions, haplotypes) :
    """for each region in phased bed, find the overlapping haplotype block"""
    regions_to_phase = []

    for r in regions : # for each maximally phased region
        min_length = r[2]-r[1] # minimum length to reach for considering an overlap

        for h in haplotypes :
            if r[0] != h[0] :
                continue
            else :
                length = get_overlap(r[1], r[2], h[1], h[2])
                if length == min_length : # very important for no tails and heads in alignments
                    regions_to_phase.append(r)
                else :
                    continue
    return regions_to_phase

def make_regions_df_for_merge(regions_to_phase, lengths) :
    """For each region to phase (sample specific), make a dataframe with coordinates and whether it is an haplotype;"""
    dfs = []

    for chrom, lgt in lengths.items() :
        chrs = np.array([chrom for i in range(lgt)])
        cpos = np.arange(1,lgt+1) # VCF is 1-indexed
        ishap = np.zeros(lgt, dtype=np.uint)
        blockn = np.zeros(lgt, dtype=np.uint)

        for i, r in enumerate(regions_to_phase) :
            if r[0] == chrom :
                ishap[r[1]:r[2]] = 1
                blockn[r[1]:r[2]] = i

        dfs.append(pd.DataFrame().from_dict({"chr":chrs, "pos":cpos, "ishap":ishap, "nid":blockn}))

    all_r_dfs = pd.concat(dfs)

    return all_r_dfs

##################### WRITING DATA #####################

def write_new_fasta(output, seq, only_phased) :
    """From a set of haplotype blocks, outputs a fasta file"""
    f = open(output, "w")

    new_seqs = {}
    for n, blocks in only_phased.groupby("nid") :
        start = blocks["pos"].min()
        end = blocks["pos"].max()
        ref = blocks["chr"].iloc[0]
        # start and end here are 1-based but python is 0-based and end is excluded
        # So we take row["start"]-1 to be exactly on the region (in VCF 1 == 0 in python so ~ 0 == -1)
        region = seq[ref][start-1:end]

        # start and end here are 1-based and end is included (for VCF and GFF)
        new_name = "{}_{}_{}".format(ref, start, end)
        new_seqs[new_name] = region

        f.write(">{}\n{}\n".format(new_name, region))

    f.close()

    return new_seqs

def write_new_vcf(output, phased_vcf, new_sequences, only_phased) : # Need pysam for this function (renaming records in VCF)
    """Writes a .VCF file containing only phasing-relevant information"""
    vcf_in = VariantFile(phased_vcf)  # auto-detect input format
    vcf_out = VariantFile(output, 'w', header=vcf_in.header, threads=1)

    # Remove contigs not in vcf anymore
    previous_contigs = list(vcf_in.header.contigs.keys())
    for ctg in previous_contigs :
        vcf_out.header.contigs.remove_header(ctg)

    # Add new contigs to vcf
    for name, seq in new_sequences.items() :
        vcf_out.header.contigs.add(name, length=len(seq))

    dc = {"ctg":[], "pos":[], "record":[]}
    for record in vcf_in :
        dc["ctg"].append(record.chrom)
        dc["pos"].append(record.pos)
        dc["record"].append(record)

    vcf_in.close()

    vcf_df = pd.DataFrame().from_dict(dc) #, dtype={"ctg":"category", "pos":"int"})
    del dc # frees memory

    new_blocks = [] # list of pandas dataframe

    gb = only_phased.groupby("nid")

    total = len(gb)
    tenpercent = int(0.1*total)
    i = 0
    for n, blocks in gb :
        if i % tenpercent == 0 :
            print("Processed {} / {} blocks".format(i, total))

        start = blocks["pos"].min()
        end = blocks["pos"].max()
        chrom = blocks["chr"].iloc[0]

        new_sequence_name = "{}_{}_{}".format(chrom, start, end)
        lines = vcf_df.query("ctg == @chrom & pos >= @start & pos <= @end")
        # update the coordinates
        data = lines.copy()
        data.loc[:]["pos"] = data.loc[:]["pos"] - (start-1) # since start is 1-based, we must substract start-1
        data.loc[:]["ctg"].replace({chrom:new_sequence_name}, inplace=True)
        new_blocks.append(data)

        i += 1

    log("Writing new .vcf file...")
    for i, block in enumerate(new_blocks) :
        if i % tenpercent == 0 :
            print("Processed {} / {} blocks".format(i, total))
        #print(block)
        for n, row in block.iterrows() :
            #print(row)
            nr = vcf_out.new_record()
            nr = row["record"].copy()
            nr.translate(vcf_out.header)
            nr.pos = row["pos"]
            nr.chrom = row["ctg"]
            vcf_out.write(nr)

    vcf_out.close()

def run_bcftools(output, ref, vcf) :
    # 1. convert the vcf to bcf
    bcf = output + ".bcf"
    cmd = "bcftools convert -o {} -O b {}".format(bcf, vcf)
    print("\tConverting: {}".format(cmd))
    run(cmd)

    # 2. index the newly created .bcf
    cmd = "bcftools index -f {}".format(bcf)
    print("\tIndexing: {}".format(cmd))
    run(cmd)

    # 3. sort the bcf
    sorted_bcf = output + ".sorted.bcf"
    cmd = "bcftools sort -O b -o {} {}".format(sorted_bcf, bcf)
    print("\tSorting: {}".format(cmd))
    run(cmd)

    # 4. index the vcf with updated coordinates
    cmd = "bcftools index -f {}".format(sorted_bcf)
    print("\tIndexing: {}".format(cmd))
    run(cmd)

    # 5. Output the haplotypes
    hap1_fasta_out = "{}.hap1.fa".format(output)
    cmd = "bcftools consensus -H 1 -o {} -f {} {}".format(hap1_fasta_out, ref, sorted_bcf)
    print("\tPhasing: {}".format(cmd))
    run(cmd)

    hap2_fasta_out = "{}.hap2.fa".format(output)
    cmd = "bcftools consensus -H 2 -o {} -f {} {}".format(hap2_fasta_out, ref, sorted_bcf)
    print("\tPhasing: {}".format(cmd))
    run(cmd)

##################### MAIN #####################

def phase(args) :
    """Parse args and run program"""

    # First check that bcftools is in $PATH
    if not which("bcftools") :
        print("Missing dependency: bcftools")
        sys.exit(1)

    # IO
    reference = args.FASTA[0]
    bed_file = args.BED[0]
    phase_outdir = args.VCF[0]
    outdir = args.outdir[0]

    if not os.path.exists(outdir) :
        os.makedirs(outdir) # Create directory following path
    else :
        log("WARNING: Output directory already exists: {}.".format(outdir))

    print("# Phasing tool")
    print("Haploid reference:\t{0}".format(reference))
    print("Regions to phase:\t{0}".format(bed_file))
    print("Phase vcf output directory:\t{0}".format(phase_outdir))
    print("Output directory path:\t\t{0}".format(outdir))

    log("Reading reference...") # common
    refseqs = parse_fasta_generic(reference)
    lengths = {k:len(v) for k,v in refseqs.items()}

    log("Parsing .bed file...") # common
    regions = parse_bed(bed_file)

    log("Parsing samples...") # common
    samples = parse_samples(phase_outdir) # format: {sample: {vcf:vcf, blocks:blocks}}

    i = 1
    t = len(samples)
    for sm, input_files in samples.items() :

        phased_blocks_hapcut2 = input_files["blocks"]
        phased_vcf = input_files["vcf"]

        log("{} ({}/{}): Parsing phased blocks file...".format(sm, i, t)) # sample wise
        df, haplotypes = parse_haplotype_blocks(phased_blocks_hapcut2)

        log("{} ({}/{}): Checking phaseable regions...".format(sm, i, t)) # sample wise
        regions_to_phase = find_overlaps(regions, haplotypes)

        log("{} ({}/{}): Getting regions to phase...".format(sm, i, t)) # sample wise
        all_r_df = make_regions_df_for_merge(regions_to_phase, lengths)
        all_r_df_only_hap = all_r_df.query("ishap > 0")
        only_phased = pd.merge(df, all_r_df_only_hap, on=["chr", "pos"], how="right")

        log("{} ({}/{}): Outputting haploid fasta of regions to phase...".format(sm, i, t)) # sample wise
        #outfa = "{}.fa".format(output)
        fa_out = os.path.join(outdir, sm+".fa")
        new_seqs = write_new_fasta(fa_out, refseqs, only_phased)

        log("{} ({}/{}): Writing a new phased .VCF with per-region information...".format(sm, i, t)) # sample wise
        #vcf_out = "{}.vcf".format(output) # Need be uncompressed
        vcf_out = os.path.join(outdir, sm + ".phased.vcf") # Need be uncompressed
        write_new_vcf(vcf_out, phased_vcf, new_seqs, only_phased)

        log("{} ({}/{}): Running bcftools consensus...".format(sm, i, t)) # sample wise
        run_bcftools(output, fa_out, vcf_out)

        i += 1
