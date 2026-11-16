#!/usr/bin/env python3
import sys
from syn_nonsyn_from_gff import (read_fasta, parse_gff3_cds,
    choose_longest_tx_per_gene, build_cds_seq_and_posmap, CODON_TABLE)

gff, fasta, kappa, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]

genome = read_fasta(fasta)
cds_by_tx, gene_by_tx, strand_by_tx = parse_gff3_cds(gff)
best = choose_longest_tx_per_gene(cds_by_tx, gene_by_tx)

transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

with open(out, "w") as o:
    o.write("gene\ttranscript\tnonsyn_w\tsyn_w\n")
    for gene, tx in best.items():
        blocks = cds_by_tx[tx]
        if len(set(b[0] for b in blocks)) != 1:
            continue
        chrom = blocks[0][0]
        if chrom not in genome:
            continue
        seq, _ = build_cds_seq_and_posmap(genome[chrom], blocks, strand_by_tx[tx])
        nonsyn = 0.0
        syn = 0.0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
            for j in range(3):
                for b in "ACGT":
                    if b == codon[j]:
                        continue
                    new = codon[:j] + b + codon[j + 1:]
                    aa2 = CODON_TABLE.get(new)
                    if aa2 is None or aa2 == "*":
                        continue
                    w = kappa if (codon[j], b) in transitions else 1.0
                    if aa2 == aa:
                        syn += w
                    else:
                        nonsyn += w
        o.write(f"{gene}\t{tx}\t{nonsyn}\t{syn}\n")
