"""
Converts autoPVS1 results to vcf file.
Done so that an existing tool can be used to augment an existing vcf file
"""

import argparse
import subprocess
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--autopvs1_tsv', help='autoPVS1 file to convert')
parser.add_argument('--output_basename', help='output_basename')
parser.add_argument('--header', action='store_true', help='Add if tsv file has header')

args = parser.parse_args()

# General input looks like this:
# vcf_id	SYMBOL	Feature	trans_name	consequence	strength_raw	strength	criterion
# 1-925952-G-A	SAMD11	ENST00000616016	not_lof	missense_variant	Unmet	Unmet	na
# 1-925956-C-T	SAMD11	ENST00000616016	not_lof	synonymous_variant	Unmet	Unmet	na
# 1-925969-C-T	SAMD11	ENST00000616016	not_lof	missense_variant	Unmet	Unmet	na

vcf_head = """##fileformat=VCFv4.2
##INFO=<ID=AUTOPVS1_STRENGTH_RAW,Number=.,Type=String,Description="autoPVS1 strength descriptor">
##INFO=<ID=AUTOPVS1_CRITERION,Number=.,Type=String,Description="autoPVS1 label code">
#CHROM\tPOS\tID	REF\tALT\tQUAL\tFILTER\tINFO
"""
in_file = args.autopvs1_tsv
out_file = args.output_basename + ".vcf.gz"
# open bgzip pipe for compression
p = subprocess.Popen("bgzip -c > " + out_file, shell=True, stdin=subprocess.PIPE)
p.stdin.write(bytes(vcf_head, 'utf-8'))
tsv = open(in_file)
if args.header:
    skip = next(tsv)
for line in tsv:
    info = line.rstrip('\n').split("\t")
    var_info = info[0].split("-")
    var_info.insert(2, ".")
    record = "chr" + "\t".join(var_info)
    record += "\t.\t.\tAUTOPVS1_STRENGTH_RAW=" + info[-3] + ";AUTOPVS1_CRITERION=" + info[-1] + "\n"
    p.stdin.write(bytes(record, 'utf-8'))
p.communicate()

subprocess.call("tabix " + out_file, shell=True)