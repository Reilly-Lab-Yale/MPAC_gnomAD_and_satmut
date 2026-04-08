zcat /gpfs/gibbs/pi/reilly/VariantEffects/scripts/noon_data/output_indel_consolidated/*.sorted.vcf.gz | grep -v "^#" | wc -l > total_lines.txt
