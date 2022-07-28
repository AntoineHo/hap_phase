#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Haplotype Phasing pipeline using python"""

#import errno
#import os
#import shutil
#import datetime

import argparse
import sys

def get_cov(args) :
    """Runs step 0 (get coverage per sample)"""
    from coverage import get_coverage
    get_coverage(args)

def phase_vcf(args) :
    """Runs step 1 (phase a VCF per sample)"""
    from phase_vcf import phase_vcf
    phase_vcf(args)

def filter(args) :
    """Runs step 1 (filter)"""
    from filter import filter
    filter(args)

def extract(args) :
    """Runs step 2 (extract)"""
    from extract import extract
    extract(args)

def phase(args) :
    """Runs step 3 (phase)"""
    from phase import phase
    phase(args)

def group(args) :
    """Runs step 4 (group)"""
    from group import group
    group(args)

def align(args) :
    """Runs step 5 (align)"""
    from align import align
    align(args)

def concat(args) :
    """Runs step 6 (concatenate)"""
    from concat import concat
    concat(args)

def main() :
    """Argument parser"""
    parser = argparse.ArgumentParser(description='Runs steps (filter variants, extract regions, phase haplotypes, align haplotypes and concatenate haplotypes).')
    subparsers = parser.add_subparsers(required=True, dest="filter || extract || phase || group || align || concatenate")

    # Run hapcut2
    phv = subparsers.add_parser("phase_vcf", help='Run hapcut2 using samples .bam and .vcf file')
    phv.add_argument('VCF',                  nargs=1, type=str,   help="<STRING> A phased VCF output of HapCut2.")
    phv.add_argument('BED',                  nargs=1, type=str,   help="<STRING> A path to a .bed file with regions to phase.")
    phv.add_argument('SAMPLES',              nargs=1, type=str,   help="<STRING> A path to a sample file list. Format one sample per line like: sample_name \t bam_file_path")
    phv.add_argument('-t','--threads',       nargs=1, type=int,   default=[4], required=False, help="<INT> Threads to use. Default: >= %(default)s")
    phv.add_argument('-o','--outdir',        nargs=1, type=str,   default=['out'], help="<STRING> path to the output directory. Default: %(default)s")
    phv.set_defaults(func=phase_vcf)

    #hap.add_argument('HAPCUT2',      nargs=1, type=str, help="<STRING> A path to HAPCUT2 executable.")
    #hap.add_argument('extractHAIRS', nargs=1, type=str, help="<STRING> A path to HAPCUT2 executable.")

    # Filter haplotype blocks
    fil = subparsers.add_parser("filter", help='Filter all valid haplotype blocks and outputs a .BED file.')
    fil.add_argument('FASTA',    nargs=1, type=str, help="<STRING> A fasta file containing the HAPLOID reference genome.")
    fil.add_argument('COV',      nargs=1, type=str, help="<STRING> A coverage file (from sambamba depth base) with all sites info from the reference.")
    fil.add_argument('VCF',      nargs=1, type=str, help="<STRING> A phased VCF output of HapCut2.")
    fil.add_argument('HAPCUT2',  nargs=1, type=str, help="<STRING> A phased block file output of HapCut2.")
    fil.add_argument('SAMPLE',   nargs=1, type=str, help="<STRING> Sample name.")
    fil.add_argument('-ml','--min-length',   nargs=1, type=int,   default=[50],   required=False,help="<INT> Minimum length of phased region. Default: >= %(default)s")
    fil.add_argument('-mc','--min-cov',      nargs=1, type=int,   default=[25],   required=False,help="<INT> Minimum average coverage in phased block. Default: >= %(default)s")
    fil.add_argument('-ma','--min-af',       nargs=1, type=float, default=[0.35], required=False,help="<INT> Minimum average allele frequency in phased block. Default: >= %(default)s")
    fil.add_argument('-mmq','--min-mis-qual',nargs=1, type=int,   default=[20],   required=False,help="<INT> Minimum mismatch quality in phased block. Default: >= %(default)s")
    fil.add_argument('-o','--output', nargs=1, type=str, default=['out'],help="<STRING> prefix for the output files. Default: %(default)s")
    fil.set_defaults(func=filter)

    # Extract maximum set of commonly phased haplotype blocks
    ext = subparsers.add_parser("extract", help='Find the maximum set of region phased in n samples.')
    ext.add_argument('DIR',     nargs=1, type=str, help="<STRING> Directory of the .bed files.")
    ext.add_argument('N',       nargs=1, type=int, help="<INT> Minimum number of samples phased to filter the regions.")
    ext.add_argument('REF',     nargs=1, type=str, help="<STRING> path to the fasta reference file.")
    ext.add_argument('-ml','--min-length',   nargs=1, type=int, default=[100],   help="<INT> Minimum length of phased region. Default: >= %(default)s")
    ext.add_argument('-o','--output',        nargs=1, type=str, default=['out'], help="<STRING> name for the output file. Default: %(default)s")
    ext.set_defaults(func=extract)

    # Obtain a phased vcf from haploid reference using haplotype blocks
    pha = subparsers.add_parser("phase", help='Using an haploid reference, a phased .vcf and a region .bed file, outputs phased sequences.')
    pha.add_argument('FASTA',    nargs=1, type=str, help="<STRING> A fasta file containing the HAPLOID reference genome.")
    pha.add_argument('BED',      nargs=1, type=str, help="<STRING> A bed file containing the regions to try to phase.")
    pha.add_argument('VCF',      nargs=1, type=str, help="<STRING> A phased VCF output of HapCut2.")
    pha.add_argument('HAPCUT2',  nargs=1, type=str, help="<STRING> A phased block file output of HapCut2.")
    pha.add_argument('-o','--output', nargs=1, type=str, default=['out'],help="<STRING> prefix for the output files. Default: %(default)s")
    pha.set_defaults(func=phase)

    # Get one fasta file per phased block
    gro = subparsers.add_parser("group", help='Takes a bed file and samples names and haplotypes then outputs one fasta per phased region.')
    gro.add_argument('FASTA',    nargs=1, type=str, help="<STRING> A fasta file containing the HAPLOID reference genome.")
    gro.add_argument('BED',      nargs=1, type=str, help="<STRING> A bed file containing the phased regions.")
    gro.add_argument('SAMPLES',  nargs=1, type=str, help="<STRING> Samples to compare.")
    gro.add_argument('SAMPLEDIR',nargs=1, type=str, help="<STRING> Path to directory containing samples.")
    gro.add_argument('OUTPUT',   nargs=1, type=str, help="<STRING> A path to a directory for outputting files.")
    gro.add_argument('-l','--length', nargs=1, type=int, default=[100], required=False, help="<INT> Minimum length. Default: >= %(default)s")
    gro.add_argument('-m','--min',    nargs=1, type=int, default=[20],  required=False, help="<INT> Minimum number of samples sharing the phased region. Default: >= %(default)s")
    gro.set_defaults(func=group)

    # Align using mafft
    aln = subparsers.add_parser("align", help='Take phased fasta and align them using mafft.')
    aln.add_argument('DIR',    nargs=1, type=str, help="<STRING> A path to a directory containing phased regions folders.")
    aln.add_argument('OUTPUT', nargs=1, type=str, help="<STRING> A path to a directory for outputting files.")
    aln.set_defaults(func=align)

    # Concatenate the sequences
    con = subparsers.add_parser("concat", help='Fetch all .mafft files and concatenate the sequences.')
    con.add_argument('DIR',     nargs=1, type=str, help="<STRING> Directory to search for all the .mafft files.")
    con.add_argument('-o','--output', nargs=1, type=str, default=['out'],help="<STRING> name for the output file. Default: %(default)s")
    con.add_argument('--samples', nargs='*', type=str, help="Sample names to concatenate.")
    con.set_defaults(func=concat)

    args = parser.parse_args()
    args.func(args)
    sys.exit(0)

if __name__ == '__main__':
	main()
