"""
Converts autoPVS1 results to vcf file.
Done so that an existing tool can be used to augment an existing vcf file
"""

import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--autopvs1_tsv', help='autoPVS1 file to convert')
parser.add_argument('--header', action='store_true', help='Add if tsv file has header')

args = parser.parse_args()

# General input looks like this:
# vcf_id	SYMBOL	Feature	trans_name	consequence	strength_raw	strength	criterion
# 1-925952-G-A	SAMD11	ENST00000616016	not_lof	missense_variant	Unmet	Unmet	na
# 1-925956-C-T	SAMD11	ENST00000616016	not_lof	synonymous_variant	Unmet	Unmet	na
# 1-925969-C-T	SAMD11	ENST00000616016	not_lof	missense_variant	Unmet	Unmet	na

vcf_head = """
##fileformat=VCFv4.2
##INFO=<ID=AUTOPVS1_STRENGTH_RAW,Number=.,Type=String,Description="autoPVS1 strength descriptor">
##INFO=<ID=AUTOPVS1_CRITERION,Number=.,Type=String,Description="autoPVS1 label code">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""
in_file = args.autopvs1_tsv
out_file = ".".join(infile.split('.')[:-2] + ".vcf.gz")
# open bgzip pipe for compression
p = subprocess.Popen("bgzip -c > " + out_file, shell=True, stdin=subprocess.PIPE)
p.stdin.write(vcf_head)
tsv = open(args.autopvs1_tsv)
if args.header:
    skip = next(tsv)
for line in tsv:
    info = line.rstrip('\n').split("\t")
    record = "chr" + info[0].replace("-","\t")
    record += "\t.\t.\tAUTOPVS1_STRENGTH_RAW=" + info[-3] + ";AUTOPVS1_CRITERION=" + info[-1] + "\n"
    p.stdin.write(record)
p.communicate()
