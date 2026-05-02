#!/bin/bash
set -e  # stop if anything fails

# Activate environment
module load miniconda3
conda activate /hb/groups/bernardi_lab/programs/vcf2phylip

# Convert VCF to FASTA
python /hb/groups/bernardi_lab/programs/vcf2phylipdir/vcf2phylip.py \
  -i GMIall_filtered.vcf.gz \
  --output-folder . \
  --fasta

# Load IQ-TREE
module load iq-tree/2.2.2.6

# Run tree inference
iqtree2 \
  -s GMIall_filtered.min4.fasta \
  --seqtype DNA \
  -m GTR+G4 \
  -B 1000 \
  -T AUTO