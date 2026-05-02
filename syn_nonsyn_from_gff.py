#!/usr/bin/env python3
import argparse
from collections import defaultdict

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

def translate_codon(codon: str) -> str:
    codon = codon.upper()
    return CODON_TABLE.get(codon, "X")

def read_fasta(fa_path: str) -> dict:
    seqs = {}
    name = None
    buf = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs

def parse_gff3_cds(gff_path: str):
    # cds_by_tx[tx] = list of (chrom, start1, end1, strand)
    cds_by_tx = defaultdict(list)
    gene_by_tx = {}
    strand_by_tx = {}
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = parts
            if feature != "CDS":
                continue
            start_i = int(start); end_i = int(end)

            attr_map = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_map[k] = v

            tx = attr_map.get("Parent")
            gene = attr_map.get("gene") or attr_map.get("Name") or attr_map.get("ID")
            if not tx:
                continue

            cds_by_tx[tx].append((chrom, start_i, end_i, strand))
            gene_by_tx[tx] = gene
            strand_by_tx[tx] = strand
    return cds_by_tx, gene_by_tx, strand_by_tx

def build_cds_seq_and_posmap(genome_seq: str, blocks, strand: str):
    # Sort blocks in transcription order
    if strand == "+":
        blocks_sorted = sorted(blocks, key=lambda x: x[1])
    else:
        blocks_sorted = sorted(blocks, key=lambda x: x[1], reverse=True)

    cds_parts = []
    posmap = {}  # genomic 1-based -> cds 0-based
    cds_i = 0

    for (_chrom, s, e, _st) in blocks_sorted:
        frag = genome_seq[s-1:e]
        if strand == "-":
            frag_rc = revcomp(frag)
            # map genomic positions to CDS indices after reverse-complement
            for p in range(s, e+1):
                posmap[p] = cds_i + (e - p)
            frag = frag_rc
        else:
            for p in range(s, e+1):
                posmap[p] = cds_i + (p - s)

        cds_parts.append(frag)
        cds_i += (e - s + 1)

    cds_seq = "".join(cds_parts).upper()
    return cds_seq, posmap

def read_fixed_alleles(path: str, keep_ids: set):
    recs = defaultdict(list)  # contig -> list of dict
    with open(path) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            chrom, pos, ref, alt, pac, gulf = line.strip().split()[:6]
            pos = int(pos)
            sid = f"{chrom}:{pos}"
            if sid not in keep_ids:
                continue
            recs[chrom].append({
                "chrom": chrom, "pos": pos, "sid": sid,
                "ref": ref.upper(), "alt": alt.upper(),
                "pac": pac.upper(), "gulf": gulf.upper()
            })
    for c in recs:
        recs[c].sort(key=lambda r: r["pos"])
    return recs

def choose_longest_tx_per_gene(cds_by_tx, gene_by_tx):
    best = {}
    best_len = {}
    for tx, blocks in cds_by_tx.items():
        gene = gene_by_tx.get(tx)
        if not gene:
            continue
        L = sum((e - s + 1) for (_c, s, e, _st) in blocks)
        if gene not in best_len or L > best_len[gene]:
            best_len[gene] = L
            best[gene] = tx
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fixed", required=True)   # PACvsGULF.fixed_alleles.tsv
    ap.add_argument("--ids", required=True)     # fixed_in_CDS.ids.txt
    ap.add_argument("--out", required=True)     # output TSV
    args = ap.parse_args()

    keep_ids = set(line.strip() for line in open(args.ids) if line.strip())

    genome = read_fasta(args.fasta)
    cds_by_tx, gene_by_tx, strand_by_tx = parse_gff3_cds(args.gff)
    fixed_by_contig = read_fixed_alleles(args.fixed, keep_ids)

    # pick one transcript per gene (longest CDS)
    best_tx_for_gene = choose_longest_tx_per_gene(cds_by_tx, gene_by_tx)

    # index transcripts by contig (and require single-contig transcripts)
    tx_on_contig = defaultdict(list)
    for gene, tx in best_tx_for_gene.items():
        chroms = sorted(set(b[0] for b in cds_by_tx[tx]))
        if len(chroms) != 1:
            continue
        chrom = chroms[0]
        tx_on_contig[chrom].append((gene, tx))

    with open(args.out, "w") as out:
        out.write("\t".join([
            "sid","chrom","pos","gene","transcript","strand",
            "codon_index","codon_pos_in_codon",
            "pac_codon","gulf_codon","pac_aa","gulf_aa","class"
        ]) + "\n")

        # For each contig, build CDS+mapping per transcript and classify its SNPs
        for chrom, snps in fixed_by_contig.items():
            if chrom not in genome or chrom not in tx_on_contig:
                continue

            for gene, tx in tx_on_contig[chrom]:
                blocks = cds_by_tx[tx]
                strand = strand_by_tx.get(tx, "+")
                cds_seq, posmap = build_cds_seq_and_posmap(genome[chrom], blocks, strand)

                for r in snps:
                    pos = r["pos"]
                    if pos not in posmap:
                        continue

                    idx = posmap[pos]
                    codon_start = (idx // 3) * 3
                    codon_end = codon_start + 3
                    if codon_end > len(cds_seq):
                        cls = "out_of_frame"
                        out.write("\t".join(map(str, [
                            r["sid"], chrom, pos, gene, tx, strand,
                            (idx // 3) + 1, (idx % 3) + 1,
                            ".", ".", ".", ".", cls
                        ])) + "\n")
                        continue

                    ref_codon = cds_seq[codon_start:codon_end]
                    # build PAC/GULF codons by substituting the allele at idx
                    pac_codon = list(ref_codon)
                    gulf_codon = list(ref_codon)
                    pac_codon[idx - codon_start] = r["pac"]
                    gulf_codon[idx - codon_start] = r["gulf"]
                    pac_codon = "".join(pac_codon)
                    gulf_codon = "".join(gulf_codon)

                    pac_aa = translate_codon(pac_codon)
                    gulf_aa = translate_codon(gulf_codon)

                    if "N" in pac_codon or "N" in gulf_codon or pac_aa == "X" or gulf_aa == "X":
                        cls = "ambiguous"
                    elif pac_aa == gulf_aa:
                        cls = "synonymous"
                    else:
                        cls = "nonsynonymous"

                    out.write("\t".join(map(str, [
                        r["sid"], chrom, pos, gene, tx, strand,
                        (idx // 3) + 1, (idx % 3) + 1,
                        pac_codon, gulf_codon, pac_aa, gulf_aa, cls
                    ])) + "\n")

if __name__ == "__main__":
    main()

