#!/usr/bin/env bash
set -euo pipefail

MAF="${1:?Usage: $0 <file.maf>}"
THREADS=8

echo "Calling variants on MAF file: $MAF"
export MAF

mkdir -p vcfs

# --- Helper: inject ##contig lines into VCF header so bcftools can parse/sort/index ---
add_contigs_to_header () {
  local in_vcf="$1"
  local out_vcf="$2"

  # Collect contig names from CHROM column (skip headers)
  awk '!/^#/ {print $1}' "$in_vcf" | sort -u > "${in_vcf}.contigs"

  # Write header, insert contig lines before #CHROM, then write body
  awk -v contigs_file="${in_vcf}.contigs" '
    BEGIN {
      while ((getline c < contigs_file) > 0) contigs[++n]=c
      close(contigs_file)
    }
    /^#CHROM/ {
      for (i=1; i<=n; i++) print "##contig=<ID=" contigs[i] ">"
      print
      next
    }
    { print }
  ' "$in_vcf" > "$out_vcf"

  rm -f "${in_vcf}.contigs"
}
export -f add_contigs_to_header  

# Get sample names from MAF 
mapfile -t SAMPLES < <(grep "^s" "$MAF" | cut -f2 | cut -f1 -d'#' | sort -u)

parallel --jobs "$THREADS" --halt now,fail=1 '
  set -euo pipefail

  SAMPLE={};

  raw="vcfs/${SAMPLE}.raw.vcf"
  fixed="vcfs/${SAMPLE}.vcf"
  gz="vcfs/${SAMPLE}.vcf.gz"

  wgatools call \
    --query-regex "^${SAMPLE}#.*" \
    --sample "$SAMPLE" \
    --snp \
    --svlen 0 \
    -r "$MAF" \
    -o "$raw"

  add_contigs_to_header "$raw" "$fixed"
  rm -f "$raw"

  bcftools sort -Oz -o "$gz" "$fixed"
  bcftools index -f "$gz"
  rm -f "$fixed"
' ::: "${SAMPLES[@]}"

bcftools merge vcfs/*.vcf.gz -o merged.vcf.gz -O z --threads "$THREADS"
bcftools index -f merged.vcf.gz

bcftools annotate --remove FORMAT/QI merged.vcf.gz -o merged.vcf
sed -E -i.bak 's#\./\.#0#g; s#([0-9])\|\1#1#g' merged.vcf

rm -rf vcfs
rm -f merged.vcf.gz merged.vcf.gz.csi merged.vcf.bak



