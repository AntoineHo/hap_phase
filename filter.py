#!/usr/bin/python
# -*- coding: utf-8 -*-

# Standard library
import os
import sys
import subprocess

# Packages
import pandas as pd
import numpy as np
from pysam import VariantFile

# Modules
from utils import log
from parsers import parse_fai, parse_bed, parse_samples

# We have various informations from different files for !one! sample:
# - Reference genome
# - All base coverage from sambamba depth base
# - Hapcut2 haplotype blocks
# - Hapcut2 VCF file

# We want to obtain the coordinates of the phased regions after various filters
# - average AF in the region (VCF data)
# - coverage in the region (sambamba depth base)
# - length of the region
# - mismatch score (hapcut2)
# For instance; blocks output have parameters such as: af >= 0.3 and af <= 0.5 / coverage >= 25x / length >= 50bp / mismatch quality >= 20

# For one given sample, there will be one filtered .bed file output

##################### PARSING FILES #####################
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
    return df

def parse_vcf(input, sample) :
    """Parse the VCF file"""
    vcf_in = VariantFile(input)  # auto-detect input format
    dc = {"chr":[], "pos":[], "GT":[], "DP":[], "AF":[]}

    for i, record in enumerate(vcf_in) :
        gt = record.samples[sample]["GT"]
        dp = record.samples[sample]["DP"]
        ad = record.samples[sample]["AD"]

        if gt[0] != gt[1] :
            #print(gt)
            dc["chr"].append(record.chrom)
            dc["pos"].append(record.pos)
            dc["GT"].append("{}|{}".format(*gt))
            dc["DP"].append(dp)
            try :
                dc["AF"].append(min(ad)/sum(ad))
            except :
                dc["AF"].append(None)

    vcf_df = pd.DataFrame().from_dict(dc) #, dtype={"ctg":"category", "pos":"int"})
    return vcf_df

def parse_coverage(input) :
    cov = pd.read_csv(
        input, usecols=range(3), names=["chr", "pos", "cov"], sep="\t",
        dtype={"chr":"category", "pos":"int", "cov":"int"}, header=0, compression='infer',
    )

    cov["pos"] += 1 # change to a 1-based coordinate system
    return cov

##################### DATA TRANSFORMATION #####################

def merge_information(blocks, vcf, cov, lengths) :
    """Merge coverage file with block positions"""

    # Merge blocks and vcf data
    blocks_vcf = pd.merge(blocks, vcf, on = ["chr", "pos"], how="inner")

    # Groupby bid and get coordinates
    block_coords = blocks_vcf.groupby(by="bid").agg({"chr":"first", "pos":["min","max"]})
    block_coords.columns = ["_".join(c for c in col) for col in block_coords.columns]
    block_coords = block_coords.reset_index()

    # Create vectors with positions and is a bid or not
    in_blocks = {chrom:np.zeros(lg, dtype=int) for chrom, lg in lengths.items()}
    positions = {chrom:np.arange(1, lg+1, 1, dtype=int) for chrom, lg in lengths.items()}
    chrs = {chrom:[chrom for i in range(lg)] for chrom, lg in lengths.items()}
    for n, row in block_coords.iterrows() :
        bid = row["bid"]
        chrom = row["chr_first"]
        start = row["pos_min"]
        end = row["pos_max"]
        in_blocks[chrom][start:end+1] = bid

    # Concatenate
    positions_and_bid = pd.DataFrame({"chr":[], "pos":[], "bid":[]})
    for chrom, bid in in_blocks.items() :
        pos = positions[chrom]
        chrls = chrs[chrom]
        positions_and_bid = pd.concat([positions_and_bid, pd.DataFrame({"chr":chrls, "pos":pos, "bid":bid})])

    # Merge blocks and coverage data
    positions_bid_cov = pd.merge(positions_and_bid, cov, on=["chr", "pos"], how="inner")
    # Get mean coverage in each block
    mean_cov_in_bid = positions_bid_cov.groupby(by="bid").agg({"bid":"first","cov":"mean"})
    mean_cov_in_bid = mean_cov_in_bid.reset_index(drop=True)

    # Merge coverage information and vcf information
    final_dataset = pd.merge(blocks_vcf, mean_cov_in_bid, on=["bid"], how="inner")

    return final_dataset

def filter_information(df, min_af, min_cov, min_mis, min_len) :

    group_data = df.groupby(by="bid").agg({"bid":"first", "chr":"first", "pos":["min", "max"], "AF":"mean", "cov":"mean", "mis":"min", "inforeads":"mean", "al1":"min", "al2":"min"})
    group_data.columns = ["_".join(c for c in col) for col in group_data.columns]
    group_data = group_data.reset_index(drop=True)
    group_data = group_data.assign(length=group_data.apply(lambda x: x["pos_max"]-x["pos_min"], axis="columns"))

    valid_blocks = group_data.query("AF_mean >= @min_af & cov_mean >= @min_cov & al1_min >= 0 & al2_min >= 0 & mis_min >= @min_mis & length >= @min_len")
    return valid_blocks

##################### WRITE BED #####################

def write_bed(valid, output) :
    valid[["chr_first", "pos_min", "pos_max"]].to_csv(output, sep="\t", index=False, header=False)

##################### MAIN #####################

def filter(args) :
    # IO
    reference = args.FASTA[0]
    sample_file = args.SAMPLES[0]
    coverage_directory = args.COVDIR[0]
    hapcut2_directory = args.PHASEDIR[0]
    outdir = args.outdir[0]
    #coverage_file = args.COV[0]
    #phased_blocks_hapcut2 = args.HAPCUT2[0]
    #phased_vcf = args.VCF[0]

    # PARAMS
    min_coverage = args.min_cov[0]
    min_af = args.min_af[0]
    min_mis_qual = args.min_mis_qual[0]
    min_len = args.min_length[0]

    if not os.path.exists(outdir) :
        os.makedirs(outdir) # Create directory following path
    else :
        log("WARNING: Output directory already exists: {}.".format(outdir))

    print("# Phasing tool")
    print("Haploid reference:\t{0}".format(reference))
    print("Input sample list:\t{0}".format(sample_file))
    print("Coverage directory:\t{0}".format(coverage_directory))
    print("Hapcut2 phasing directory:\t{0}".format(hapcut2_directory))
    print("Output directory path:\t\t{0}".format(output))
    print("Other parameters:")
    print("\t- Minimum coverage: {}".format(min_coverage))
    print("\t- Minimum AF: {}".format(min_af))
    print("\t- Minimum mismatch quality: {}".format(min_mis_qual))
    print(
        "===============================================================================\n"
    )

    #print("Coverage of assembly:\t{0}".format(coverage_file))
    #print("HapCut2 phased VCF:\t{0}".format(phased_vcf))
    #print("HapCut2 phased blocks:\t{0}".format(phased_blocks_hapcut2))

    log("Parsing reference index...")
    try :
        lengths = parse_fai(reference + ".fai")
    except :
        raise Exception("Could not find fasta index (.fai) file.")

    samples = parse_samples(sample_file)
    i = 1
    for sample, bam in samples.items() :

        log("Sample {}/{}: {}".format(i, len(samples), sample))

        # Fetch files
        coverage_file = os.path.join(coverage_directory, bam + ".cov.gz")
        phased_vcf = os.path.join(hapcut2_directory, sample + ".hapcut2.phased.VCF")
        block_files = os.path.join(hapcut2_directory, sample + ".hapcut2")
        output = os.path.join(outdir, sample + ".phasing.bed")

        log("Parsing coverage...")
        cov = parse_coverage(coverage_file)

        log("Parsing VCF...")
        vcf = parse_vcf(phased_vcf, sample)

        log("Parsing haplotype blocks...")
        blocks = parse_haplotype_blocks(block_files)

        log("Merging data...")
        data = merge_information(blocks, vcf, cov, lengths)

        log("Filtering data...")
        valid = filter_information(data, min_af, min_coverage, min_mis_qual, min_len)

        log("Outputting valid blocks coordinates...")
        write_bed(valid, output)

        log("Done!")

        i += 1
