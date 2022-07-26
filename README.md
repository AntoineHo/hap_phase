# hap_phase
Phase haplotypes from at least n samples using HapCut2 and a VCF

### Usage

```text
python main.py --help
usage: main.py [-h] {filter,extract,phase,group,align,concat} ...

Runs steps (filter variants, extract regions, phase haplotypes, align and concatenate them).

positional arguments:
  {filter,extract,phase,group,align,concat}
    filter              Filter all valid haplotype blocks and outputs a .BED file.
    extract             Find the maximum set of region phased in n samples.
    phase               Using an haploid reference, a phased .vcf and a region .bed file, outputs phased sequences.
    group               Takes a bed file and samples names and haplotypes then outputs one fasta per phased region.
    align               Take phased fasta and align them using mafft.
    concat              Fetch all .mafft files and concatenate the sequences.

optional arguments:
  -h, --help            show this help message and exit
```
