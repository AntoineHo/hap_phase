#!/usr/bin/python
# -*- coding: utf-8 -*-


def parse_bed(bed) :
    """Reads a .bed file and stores coordinate information"""
    # WARN: .bed files are 0-based!
    regions = []

    regid = 1
    for n, line in enumerate(open(bed, "r")) :
        s = line.strip().split()
        regions.append((s[0], int(s[1]), int(s[2]), regid))
        regid += 1

    return regions

def parse_fai(faifile) :
    lengths = {}
    for line in open(faifile, "r") :
        s = line.strip().split("\t")
        lengths[s[0]] = int(s[1])
    return lengths

def parse_fasta_generic(reference) :
    seq = {}

    f = open(reference, "r")
    for line in f :
        if line[0] == ">" :
            name = line[1:].strip().split()[0] # remove the end at first space
            seq[name] = ""
        else :
            seq[name] += line.strip()

    f.close()
    return seq

def parse_fasta_specific(fasta) :
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
