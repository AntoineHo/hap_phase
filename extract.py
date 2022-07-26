#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import pandas as pd

from utils import log
from parsers import parse_fai, parse_bed

# We have per-sample bed files that store all regions that are phased in a given sample
# In these .bed files the phased regions are filtered:
# af >= 0.3 <= 0.5 / coverage >= 25x / length >= 50bp / mismatch quality >= 20

# ====================== WARNING =====================
# The bed files were obtained with filter_phased_regions so the .bed                !!!
# have identical positions with the VCF file (1-based) but .bed should be 0-based   !!!
# So the output of this script must be 0-based to work with the phase_bed.py script !!!
# ====================================================

# Now the problem is to find an optimal phased set in which as many sequences in
# as many samples are phased. The idea is to loop through .bed files and add +1
# on all bases from a region if it is phased in the current bed file
# Then we can select all the regions that are phased in at least "n" samples by
# filtering this count to keep bases with a score of at least n

# Now that we have a new "max_phased.bed" file we need to loop through samples .bed
# again and take the intersection.
# We need an asymetric intersection : keeping only regions fully phased if they
# are shorter than in the "max_phased" bed file then we do not take them at all
# <- avoid large missing portions at the end of the region

# Now we need to use the intersected sample .bed file to phase the sequences

# We have to add one step here that will create Ns in the region if it is in a
# "repeats.bed" file. This will avoid counting phased heterozygous sites in the
# region like confident sites to do so we need to look through the VCF and change
# heterozygous 0/1 (whathever 0 and 1 may be) to N/N sites so that when
# bcftools consensus is used it will output N in both haplotypes

# I tested functions in .ipynb and toy datasets before writing a script that will
# do all these steps automatically

def maximally_phased_dataset(bed_dir, output, lengths, n, min_length=100) :
    """
    Reads per sample .bed file and outputs regions of at least @min_length
    phased in at least @n samples
    """
    arrays = {k:np.zeros(lg) for k, lg in lengths.items()} # 0 - indexed because python (like bed files)

    bed_files = [os.path.join(bed_dir, file) for file in os.listdir(bed_dir)] # list of filepaths from directory
    bed_files = [f for f in bed_files if f.endswith(".bed")] # remove non-bed files

    if len(bed_files) < n :
        raise Exception("Not enough files to pass the filters on regions: output will be empty.")

    for i, file in enumerate(bed_files) :
        regions = parse_bed(file)
        for r in regions :
            arrays[r[0]][r[1]:r[2]] += 1

    new_bed_data = []
    for chrom, regions in arrays.items() :
        pos = np.argwhere(regions >= n) # (0-based) positions of regions phased in more than min_sample
        transpos = pos.transpose()[0] # linearize new array

        # Groupby contiguous values
        groups = np.split(transpos, np.where(np.diff(transpos) != 1)[0]+1)
        # Get start, end and length of each region
        coords_length = np.array([get_coords_and_length(x) for x in groups])
        # Convert to pandas dataframe
        cl = pd.DataFrame(coords_length, columns=["start", "end", "length"])
        cl = cl.astype(np.uint)
        # Keep only regions longer than @min_length
        cl = cl.query("length >= @min_length")

        chrom_col = [chrom for i in range(len(cl))]
        cl = cl.assign(chrom=chrom_col)
        new_bed_data.append(cl[["chrom", "start", "end"]])

    # Concatenate all chromosomes
    new_bed = pd.concat(new_bed_data)
    # restore 0-based bed file
    new_bed["start"] -= 1
    new_bed["end"] -= 1
    # Output coordinates to .tsv
    new_bed.to_csv(output, sep="\t", header=False, index=False)
    return

def get_coords_and_length(x) :
    """Returns min, max and length of a ndarray of continuous values"""
    return np.array([np.min(x), np.max(x), np.size(x)], dtype=np.uint)

def extract(args) :

    # IO
    directory = args.DIR[0]
    reference = args.REF[0]
    min_sample = args.N[0]
    min_len = args.min_length[0]
    #repeats = args.repeats[0]
    output = args.output[0]
    #samples = args.samples

    log("Reading reference...")
    try :
        fasindex = reference + ".fai"
        open(fasindex, "r")
    except :
        raise Exception("Could not find .fai file.")
    lengths = parse_fai(fasindex)

    log("Finding maximally phased dataset in at least {} samples of at least {}bp...".format(min_sample, min_len))
    maximally_phased_dataset(directory, output, lengths, min_sample, min_len)

    log("Output per-sample regions to phase...")
    # These regions must be totally included: same size or larger than the corresponding
    # region in the maximally phased bed

    log("Done!")
