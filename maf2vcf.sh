#!/bin/bash
MAF="$1"
THREADS=8
echo "Calling variants on MAF file: $MAF"
export MAF=$MAF
mkdir -p vcfs
parallel --jobs $THREADS 'SAMPLE={}; wgatools call --query-regex "^${SAMPLE}#.*" --sample $SAMPLE --snp -o vcfs/$SAMPLE.vcf --svlen 0 -r  $MAF; bgzip vcfs/$SAMPLE.vcf; bcftools index vcfs/$SAMPLE.vcf.gz' ::: $(grep "^s" $MAF | cut -f2 | cut -f1 -d'#' | sort | uniq)
bcftools merge vcfs/*.vcf.gz -o merged.vcf.gz -O z --threads $THREADS
bcftools index merged.vcf.gz
bcftools annotate --remove FORMAT/QI  merged.vcf.gz -o merged.vcf
sed -E -i.bak 's#\./\.#0#g; s#([0-9])\|\1#1#g' merged.vcf
rm -r vcfs
rm merged.vcf.gz merged.vcf.gz.csi merged.vcf.bak

