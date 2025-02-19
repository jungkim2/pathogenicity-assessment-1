#!/usr/bin/env bash
vcf_file=$1
exp_args=("$@")

## get vcf filename
vcf_filtered_file=${vcf_file%.vcf*}."filtered.vcf"

echo "vcf file: $vcf_file ";



# default filters
cmd="bcftools filter -i 'INFO/AF >= 0.2' $vcf_file | bcftools filter -i 'INFO/DP >=15' | bcftools view -f PASS"

## loop through args for other user-defined filters
for i in "${exp_args[@]:1}"; do
  #echo filtering for... $i
  cmd+=" | bcftools filter -i '$i'"
#  cmd+=" $vcf_file "
  #echo $cmd
done

cmd+=" > $vcf_filtered_file"

echo "cmd: " $cmd
eval "$cmd"
