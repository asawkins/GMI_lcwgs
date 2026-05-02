#!/bin/bash
source config.sh


##############################
# 1) Lift over annotation to your reference
##############################
liftoff \
  "$FASTA" \
  "$REF_GENOME" \
  -g "$GFF" \
  -o liftoff.gff3 \
  -u liftoff.unmapped.txt \
  -p 16 \
  -s 0.8 \
  -a 0.8 \
  -d 1.0 \
  -flank 0.5 \
  -copies \
  -cds

##############################
# 2) Subset VCF to outlier loci
##############################
cd "$WORKDIR"

# Make BED from CSV (expects columns: contig in col2, position in col3)
awk -F',' 'NR>1 && $2!="" && $3!="" {OFS="\t"; print $2, $3-1, $3}' \
  GMIall_outlierloci_metadata.csv \
  > GMIall_outlierloci.bed

sort -k1,1 -k2,2n GMIall_outlierloci.bed > GMIall_outlierloci.sorted.bed

# Extract outlier SNPs only
bcftools view -R GMIall_outlierloci.sorted.bed -v snps -Oz \
  -o outliers.only.snps.vcf.gz \
  $VCF

bcftools index -t outliers.only.snps.vcf.gz
bcftools index -n outliers.only.snps.vcf.gz   # expect ~466

##############################
# 3) Fixed-within-pop filtering (PAC and GULF)
##############################
# Inputs: pacific_samples.txt and gulf_samples.txt
# (one sample ID per line; must match VCF header sample names)

# Pacific fixed sites (AF==0 or AF==1 within PAC)
bcftools view -S pacific_samples.txt -Ou outliers.only.snps.vcf.gz \
| bcftools +fill-tags -Ou -- -t AF,AN,AC \
| bcftools view -i 'AN>0 && (AF==0 || AF==1)' -Oz -o PAC.outliers.fixed.vcf.gz
bcftools index -t PAC.outliers.fixed.vcf.gz

# Gulf fixed sites (AF==0 or AF==1 within GULF)
bcftools view -S gulf_samples.txt -Ou outliers.only.snps.vcf.gz \
| bcftools +fill-tags -Ou -- -t AF,AN,AC \
| bcftools view -i 'AN>0 && (AF==0 || AF==1)' -Oz -o GULF.outliers.fixed.vcf.gz
bcftools index -t GULF.outliers.fixed.vcf.gz

# Split by REF-fixed vs ALT-fixed for each population
bcftools view -i 'AF==0' -Oz -o PAC.ref.vcf.gz  PAC.outliers.fixed.vcf.gz
bcftools view -i 'AF==1' -Oz -o PAC.alt.vcf.gz  PAC.outliers.fixed.vcf.gz
bcftools view -i 'AF==0' -Oz -o GULF.ref.vcf.gz GULF.outliers.fixed.vcf.gz
bcftools view -i 'AF==1' -Oz -o GULF.alt.vcf.gz GULF.outliers.fixed.vcf.gz

bcftools index -t PAC.ref.vcf.gz
bcftools index -t PAC.alt.vcf.gz
bcftools index -t GULF.ref.vcf.gz
bcftools index -t GULF.alt.vcf.gz

##############################
# 4) Fixed differences PAC vs GULF
##############################
mkdir -p isec_tmp

# sites where PAC fixed ALT and GULF fixed REF
bcftools isec -n=2 -w1 -Oz \
  -p isec_tmp/PACalt_GULFref \
  PAC.alt.vcf.gz GULF.ref.vcf.gz

# sites where PAC fixed REF and GULF fixed ALT
bcftools isec -n=2 -w1 -Oz \
  -p isec_tmp/PACref_GULFalt \
  PAC.ref.vcf.gz GULF.alt.vcf.gz

# combine both sets
bcftools concat -Oz \
  isec_tmp/PACalt_GULFref/0000.vcf.gz \
  isec_tmp/PACref_GULFalt/0000.vcf.gz \
  -o PACvsGULF.fixed_differences.vcf.gz

bcftools index -t PACvsGULF.fixed_differences.vcf.gz
bcftools index -n PACvsGULF.fixed_differences.vcf.gz   # expect 466

##############################
# 5) Build per-pop fixed allele table (PAC vs GULF)
##############################
# Extract REF/ALT/AF from each pop’s fixed set
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' PAC.outliers.fixed.vcf.gz > PAC.fixed.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' GULF.outliers.fixed.vcf.gz > GULF.fixed.tsv

# Key by chr:pos
awk 'BEGIN{OFS="\t"} {key=$1":"$2; pac_ref[key]=$3; pac_alt[key]=$4; pac_af[key]=$5}
     END{for(k in pac_af) print k, pac_ref[k], pac_alt[k], pac_af[k]}' PAC.fixed.tsv \
| sort -k1,1 > PAC.fixed.key.tsv

awk 'BEGIN{OFS="\t"} {key=$1":"$2; gulf_ref[key]=$3; gulf_alt[key]=$4; gulf_af[key]=$5}
     END{for(k in gulf_af) print k, gulf_ref[k], gulf_alt[k], gulf_af[k]}' GULF.fixed.tsv \
| sort -k1,1 > GULF.fixed.key.tsv

# Keep only sites fixed in both, but with different fixed allele (AF differs 0 vs 1)
join -1 1 -2 1 PAC.fixed.key.tsv GULF.fixed.key.tsv \
| awk 'BEGIN{OFS="\t"}
  {
    key=$1
    pac_ref=$2; pac_alt=$3; pac_af=$4
    gulf_ref=$5; gulf_alt=$6; gulf_af=$7
    if(pac_af!=gulf_af) print key,pac_ref,pac_alt,pac_af,gulf_af
  }' > PACvsGULF.fixed_sites.tsv

# Convert AF (0/1) into explicit alleles per pop + add header
awk 'BEGIN{OFS="\t"; print "CHROM","POS","REF","ALT","PAC","GULF"}
{
  split($1,a,":"); chr=a[1]; pos=a[2];
  ref=$2; alt=$3; pac_af=$4; gulf_af=$5;
  pac = (pac_af==0 ? ref : alt);
  gulf = (gulf_af==0 ? ref : alt);
  print chr,pos,ref,alt,pac,gulf
}' PACvsGULF.fixed_sites.tsv > PACvsGULF.fixed_alleles.tsv

##############################
# 6) Intersect fixed sites with CDS intervals
##############################
# BED with IDs
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2-1,$2,$1":"$2}' PACvsGULF.fixed_alleles.tsv \
  > PACvsGULF.fixed_alleles.bed

# CDS BED from liftoff GFF3 (0-based, half-open)
awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS"{print $1, $4-1, $5}' liftoff.gff3 > liftoff.CDS.bed

# Identify which fixed sites fall inside ANY CDS
awk '
  NR==FNR {cds[$1] = cds[$1] sprintf(" %d-%d", $2, $3); next}
  {
    chr=$1; pos0=$2; id=$4;
    n=split(cds[chr], a, " ");
    for(i=1;i<=n;i++){
      if(a[i]=="") continue;
      split(a[i], r, "-");
      if(pos0 >= r[1] && pos0 < r[2]) {print id; break}
    }
  }
' liftoff.CDS.bed PACvsGULF.fixed_alleles.bed > fixed_in_CDS.ids.txt

wc -l fixed_in_CDS.ids.txt   # expect 171

##############################
# 7) Synonymous vs nonsynonymous classification (Python)

# see scripts/syn_nonsyn_from_gff.py

##############################
# Save script as: syn_nonsyn_from_gff.py
chmod +x syn_nonsyn_from_gff.py

python3 syn_nonsyn_from_gff.py \
  --gff liftoff.gff3 \
  --fasta "$FASTA" \
  --fixed PACvsGULF.fixed_alleles.tsv \
  --ids fixed_in_CDS.ids.txt \
  --out fixed_in_CDS.syn_nonsyn.tsv

# Summarize counts
cut -f13 fixed_in_CDS.syn_nonsyn.tsv | tail -n +2 | sort | uniq -c


