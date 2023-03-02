#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo $script_directory
input_file="$script_directory"/""

## usage message to print for options
usage() {
  echo "Usage: $0 [-v <*.vcf>][-i <*.txt.intervar>] [-a <*autopvs1.tsv>]
  [-w <cavatica || user] [-g <gnomad_var>] [-f genomAD_AF_filter <default: 0.001] [-v variant_depth_filter default: 15 ][-r variant_AF <default: 0.2>]"
  echo ""
  echo "Options:"
  echo "  -v    vcf file"
  echo "  -i    intervar results file"
  echo "  -a    autopvs1 results file"
  echo "  -w    workflow type, must be either 'cavatica' or 'user'"
  echo "  -g    gnomAD variable, default: gnomad_3_1_1_AF_non_cancer if -w is 'cavatica'"
  echo "  -f    gnomAD allele frequency filter (default: 0.001)"
  echo "  -v    variant depth filter (default: 15)"
  echo "  -r    variant_AF <variant allele frequency (default: 0.2)"
  echo "  -c    clinvar file path (default: clinvar_20230218.vcf.gz)"
  echo "  -p    prefix for output file naming (default: 'sample')"

  echo "  -h    Display usage information."
  1>&2; exit 1; }

## default values for arguments (will be overwritten if arg passed)
clinvar_version="clinvar_20211225"
genomAD_AF_filter=0.001
variant_depth_filter=15
variant_AF=.2
prefix="sample"

## get and save arguments
while getopts ":p:v:i:a:w:g:f:v:r:c:h" arg; do
    case "$arg" in
        p) # prefix
          prefix="$OPTARG"
          ;;
        v) # vcf file
          vcf_file="$OPTARG"
          ;;
        i) # intervar file
          intervar_file="$OPTARG"
          ;;
        a) # autopvs1 file
          autopvs1_file="$OPTARG"
          ;;
        w) #workflow type
          workflow_type="$OPTARG"
          ;;
        g) #gnomad variable/column
          gnomad_var="$OPTARG"
          ;;
        f) #genomAD_AF_filter
          genomAD_AF_filter="$OPTARG"
          ;;
        v) #variant_depth_filter type
          variant_depth_filter="$OPTARG"
          ;;
        r) #variant_AF type
          variant_AF="$OPTARG"
          ;;
        c) ## clinvar version
          clinvar_version="$OPTARG"
          ;;
        h | *) # Display help.
          usage
          exit 0
          ;;
    esac
done

## check to see if workflow is either "cavatica" or "user"
if [[ -z ${workflow_type} ]]
then
  echo "ERROR: require -w option ('cavatica' or 'user')"
  exit 1;
fi

## if cavatica worklflow, save gnomad variable as "gnomad_3_1_1_AF_non_cancer" and run cmd
if [ "$workflow_type" == 'cavatica' ]
then
  gnomad_var="gnomad_3_1_1_AF_non_cancer"

  ## check if files exists and then call R script
  if [[ -f "$vcf_file" && -f "$intervar_file"  && -f "$autopvs1_file" ]];
  then
    date
    echo "Rscript 01-annotate_variants.R --vcf $vcf_file --intervar $intervar_file --autopvs1 $autopvs1_file --gnomad_variable $gnomad_var --gnomad_af $genomAD_AF_filter --variant_depth $variant_depth_filter --variant_af $variant_AF --sample_name $prefix"
    #Rscript 01-annotate_variants.R --vcf $vcf_file --intervar $intervar_file --autopvs1 $autopvs1_file --gnomad_variable $gnomad_var --gnomad_af $genomAD_AF_filter --variant_depth $variant_depth_filter --variant_af $variant_AF
  else
    echo "error: files do not exist."
    exit 1
  fi
else

  ## retrieve clinvar vcf if specified and download to input folder

  # generate full path to download
  ftp_path="ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/"$clinvar_version".vcf.gz"
  kREGEX_CLINVAR='clinvar[_/][0-9]{8}' # note use of [0-9] to avoid \d
  
  #clinvar_version="input/clinvar.vcf.gz" ## default/latest version

  ## wget clinvar file if workflow type is non-cavaita/"user" and its specified, otherwise use default clinvar db
  #if [[ $clinvar_version =~ $kREGEX_CLINVAR && $workflow_type == 'user' ]]
  #then
    ## check to see if ftp path for clinvar version exists and if so, get it
   # if [[ `wget -S --spider $ftp_path 2>&1 | grep 'Remote file exists.'` ]]; then exit_status=$?; fi
   # if [[ $exit_status == 0 ]];
   # then
    #  echo "wget -l 3 $ftp_path -P input/ wait = TRUE"
    #else
     # echo "ERROR: clinVar file $ftp_path does note exist, try again..."
      #exit 1;
   # fi
  #else
   # echo "ERROR: clinvar format error, must provide clinvar version (ie. clinvar_20211225) to download, using default"
    #exit 1;
  #fi


  ## check to see gnomad variable option entered if workflow is "user"
  if [[ "$workflow_type" == 'user' && -z ${gnomad_var} ]]
  then
    echo "ERROR: if workflow type is of type 'user', must provide -g gnomAD_var (ie. 'gnomAD_genome_ALL') ";
    exit 1;
  fi

  ## run autopvs1 and save output to input folder
  # check to see if hg38 exists, if not then download and unzip
  if [ -f "data/hg38.fa" ]
  then
    python3 ../autopvs1/autoPVS1_from_VEP_vcf.py --genome_version hg38 --vep_vcf $vcf_file  > $prefix".vcf.vep
    Rscript 02-annotate_variants_user.R --vcf $vcf_file --intervar $intervar_file --autopvs1 $autopvs1_file --clinvar $clinvar_version --gnomad_variable $gnomad_var --gnomad_af $genomAD_AF_filter --variant_depth $variant_depth_filter --variant_af variant_AF --sample_name $prefix
  else
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P data/ wait=TRUE && gunzip data/hg38.fa.gz
    python3 ../autopvs1/autoPVS1_from_VEP_vcf.py --genome_version hg38 --vep_vcf $vcf_file  > $prefix".vcf.vep
    Rscript 02-annotate_variants_user.R --vcf $vcf_file --intervar $intervar_file --autopvs1 $autopvs1_file --clinvar $clinvar_version --gnomad_variable $gnomad_var --gnomad_af $genomAD_AF_filter --variant_depth $variant_depth_filter --variant_af variant_AF --sample_name $prefix
  fi
fi
