#!/bin/bash

module load bioinfo-itrop
module load samtools/1.22.1
module load bcftools/1.22
module load htslib/1.19

set -euo pipefail

###############################################################################

# PIPELINE SNP PANEL

#

# Usage:

#   sbatch script.sh sample.bam SAMPLE_NAME

#

# Inputs:

#   reference.fa (coffea canephora genome)

#   panel_snps.txt (panel of the 28800 snp (Hamon et al., 2017)

#   sample.bam

#

###############################################################################

############################

# USER CONFIG

############################

REF="reference.fa"

PANEL="panel_snps.txt"

BAM="${1:?Usage: sbatch script.sh sample.bam SAMPLE_NAME}"

SAMPLE="${2:?Usage: sbatch script.sh sample.bam SAMPLE_NAME}"

THREADS=8

MINMAPQ=20

MINBQ=20

OUTDIR="results_${SAMPLE}"

mkdir -p "$OUTDIR"

mkdir -p "$OUTDIR/fastas_by_sample"

############################

# ENVIRONMENT

############################

echo "===== ENVIRONMENT ====="

date

hostname

command -v samtools >/dev/null || { echo "samtools missing"; exit 1; }

command -v bcftools >/dev/null || { echo "bcftools missing"; exit 1; }

command -v python3 >/dev/null || { echo "python3 missing"; exit 1; }

samtools --version | head -n 1

bcftools --version | head -n 1

############################

# CHECK INPUTS

############################

[ -f "$REF" ] || { echo "Missing $REF"; exit 1; }

[ -f "$PANEL" ] || { echo "Missing $PANEL"; exit 1; }

[ -f "$BAM" ] || { echo "Missing BAM: $BAM"; exit 1; }

############################

# INDEX REFERENCE

############################

echo "===== INDEX REF ====="

if [ ! -f "${REF}.fai" ]; then

    samtools faidx "$REF"

fi

############################

# STEP 1 - PREPARE PANEL FILES

############################

echo "===== BUILD PANEL FILES ====="

python3 << PYCODE

panel = "$PANEL"

sample = "$SAMPLE"

outdir = "$OUTDIR"

master = open(f"{outdir}/{sample}_panel_master.tsv", "w")

sites  = open(f"{outdir}/{sample}_panel_sites.tsv", "w")

bed    = open(f"{outdir}/{sample}_panel_sites.bed", "w")

vcf    = open(f"{outdir}/{sample}_panel_sites.vcf", "w")

master.write("ORDER\\tSNPID\\tCHROM\\tPOS\\tREF\\tALT\\n")

vcf.write("##fileformat=VCFv4.2\\n")

vcf.write("##source=known_panel\\n")

vcf.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n")

order = 1

with open(panel) as f:

    for raw in f:

        line = raw.strip()

        if not line:

            continue

        cols = line.split()

        if len(cols) < 4:

            continue

        snpid   = cols[0]

        alleles = cols[1]

        chrom   = cols[2]

        pos     = int(cols[3])

        if "/" not in alleles:

            continue

        ref, alt = alleles.split("/")

        master.write(f"{order}\\t{snpid}\\t{chrom}\\t{pos}\\t{ref}\\t{alt}\\n")

        sites.write(f"{chrom}\\t{pos}\\n")

        bed.write(f"{chrom}\\t{pos-1}\\t{pos}\\n")

        vcf.write(f"{chrom}\\t{pos}\\t{snpid}\\t{ref}\\t{alt}\\t.\\tPASS\\t.\\n")

        order += 1

master.close()

sites.close()

bed.close()

vcf.close()

PYCODE

bgzip -f "$OUTDIR/${SAMPLE}_panel_sites.vcf"

bcftools index -f "$OUTDIR/${SAMPLE}_panel_sites.vcf.gz"

############################

# STEP 2 - CHECK BAM

############################

echo "===== CHECK BAM ====="

if [ ! -f "${BAM}.bai" ]; then

    samtools index "$BAM"

fi

############################

# STEP 3 - GENOTYPE PANEL SNPs

############################

echo "===== BCFTOOLS MPILEUP ====="

bcftools mpileup \

    --threads "$THREADS" \

    -f "$REF" \

    -R "$OUTDIR/${SAMPLE}_panel_sites.tsv" \

    -q "$MINMAPQ" \

    -Q "$MINBQ" \

    -a FORMAT/DP,FORMAT/AD \

    "$BAM" \

| bcftools call \

    --threads "$THREADS" \

    -m \

    -Oz \

    -o "$OUTDIR/${SAMPLE}_panel.vcf.gz"

bcftools index -f "$OUTDIR/${SAMPLE}_panel.vcf.gz"

############################

# STEP 4 - EXPORT TABLES + FASTA

############################

echo "===== EXPORT RESULTS ====="

python3 << PYCODE

import gzip

import os

import re

sample = "$SAMPLE"

outdir = "$OUTDIR"

########################################

# LOAD PANEL ORDER

########################################

master = {}

order_list = []

with open(f"{outdir}/{sample}_panel_master.tsv") as f:

    next(f)

    for line in f:

        order, snpid, chrom, pos, ref, alt = line.strip().split("\\t")

        key = (chrom, pos)

        master[key] = {

            "order": int(order),

            "snpid": snpid,

            "ref": ref,

            "alt": alt

        }

        order_list.append((int(order), key))

order_list.sort()

########################################

# IUPAC TABLE

########################################

iupac = {

    frozenset(["A"]): "A",

    frozenset(["C"]): "C",

    frozenset(["G"]): "G",

    frozenset(["T"]): "T",

    frozenset(["A", "G"]): "R",

    frozenset(["C", "T"]): "Y",

    frozenset(["G", "C"]): "S",

    frozenset(["A", "T"]): "W",

    frozenset(["G", "T"]): "K",

    frozenset(["A", "C"]): "M",

    frozenset(["A", "C", "G"]): "V",

    frozenset(["A", "C", "T"]): "H",

    frozenset(["A", "G", "T"]): "D",

    frozenset(["C", "G", "T"]): "B",

    frozenset(["A", "C", "G", "T"]): "N"

}

valid = {"A", "C", "G", "T"}

########################################

# READ VCF

########################################

vcf_file = f"{outdir}/{sample}_panel.vcf.gz"

vcf_sample = None

calls = {}

with gzip.open(vcf_file, "rt") as f:

    for line in f:

        if line.startswith("##"):

            continue

        if line.startswith("#CHROM"):

            header = line.strip().split("\\t")

            vcf_sample = header[9]

            continue

        cols = line.strip().split("\\t")

        chrom = cols[0]

        pos   = cols[1]

        ref   = cols[3]

        alt   = cols[4]

        fmt   = cols[8].split(":")

        val   = cols[9]

        key = (chrom, pos)

        if key not in master:

            continue

        fi = {x: i for i, x in enumerate(fmt)}

        z = val.split(":")

        gt = z[fi["GT"]] if "GT" in fi and fi["GT"] < len(z) else "./."

        dp = z[fi["DP"]] if "DP" in fi and fi["DP"] < len(z) else "0"

        allele_map = [ref]

        if alt != ".":

            allele_map.extend(alt.split(","))

        alleles = []

        for a in re.split(r"[/|]", gt):

            if a == ".":

                continue

            try:

                idx = int(a)

                if idx < len(allele_map):

                    alleles.append(allele_map[idx].upper())

            except:

                pass

        calls[key] = {

            "GT": alleles,

            "DP": dp

        }

########################################

# GENOTYPE TABLE

########################################

with open(f"{outdir}/{sample}_ordered_genotypes.tsv", "w") as out:

    out.write("ORDER\\tSNPID\\tCHROM\\tPOS\\tREF\\tALT\\t" + sample + "\\n")

    for order, key in order_list:

        m = master[key]

        row = [

            str(order),

            m["snpid"],

            key[0],

            key[1],

            m["ref"],

            m["alt"]

        ]

        if key in calls:

            gt = calls[key]["GT"]

            row.append("/".join(gt) if len(gt) > 0 else "NA")

        else:

            row.append("NA")

        out.write("\\t".join(row) + "\\n")

########################################

# DEPTH TABLE

########################################

with open(f"{outdir}/{sample}_ordered_depth.tsv", "w") as out:

    out.write("ORDER\\tSNPID\\tCHROM\\tPOS\\t" + sample + "\\n")

    for order, key in order_list:

        m = master[key]

        row = [

            str(order),

            m["snpid"],

            key[0],

            key[1]

        ]

        if key in calls:

            row.append(calls[key]["DP"])

        else:

            row.append("0")

        out.write("\\t".join(row) + "\\n")

########################################

# BUILD FASTA

########################################

seq = []

for order, key in order_list:

    base = "N"

    if key in calls:

        gt = [x for x in calls[key]["GT"] if x in valid]

        if len(gt) > 0:

            base = iupac.get(frozenset(gt), "N")

    seq.append(base)

sequence = "".join(seq)

########################################

# FASTA GLOBAL

########################################

with open(f"{outdir}/{sample}_phylogeny_matrix.fasta", "w") as out:

    out.write(f">{sample}\\n")

    out.write(sequence + "\\n")

########################################

# FASTA SAMPLE

########################################

def clean(x):

    return re.sub(r"[^A-Za-z0-9._-]", "_", x)

fname = f"{outdir}/fastas_by_sample/{clean(sample)}.fasta"

with open(fname, "w") as out:

    out.write(f">{sample}\\n")

    out.write(sequence + "\\n")

PYCODE

############################

# DONE

############################

echo "===== DONE ====="

echo "$OUTDIR/${SAMPLE}_panel_master.tsv"

echo "$OUTDIR/${SAMPLE}_panel_sites.tsv"

echo "$OUTDIR/${SAMPLE}_panel_sites.bed"

echo "$OUTDIR/${SAMPLE}_panel_sites.vcf.gz"

echo "$OUTDIR/${SAMPLE}_panel.vcf.gz"

echo "$OUTDIR/${SAMPLE}_ordered_genotypes.tsv"

echo "$OUTDIR/${SAMPLE}_ordered_depth.tsv"

echo "$OUTDIR/${SAMPLE}_phylogeny_matrix.fasta"

echo "$OUTDIR/fastas_by_sample/${SAMPLE}.fasta"
