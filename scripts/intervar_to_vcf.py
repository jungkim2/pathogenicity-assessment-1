"""
Converts InterVar results to vcf file.
Done so that an existing tool can be used to augment an existing vcf file
"""

import argparse
import subprocess
import re
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--intervar_tsv', help='intervar file to convert')
parser.add_argument('--output_basename', help='output_basename')

args = parser.parse_args()

# General input looks like this:
#Chr	Start	End	Ref	Alt	Ref.Gene	Func.refGene	ExonicFunc.refGene	Gene.ensGene	avsnp147	AAChange.ensGene	AAChange.refGene	clinvar: Clinvar 	 InterVar: InterVar and Evidence 	Freq_gnomAD_genome_ALL	Freq_esp6500siv2_all	Freq_1000g2015aug_all	CADD_raw	CADD_phred	SIFT_score	GERP++_RS	phyloP46way_placental	dbscSNV_ADA_SCORE	dbscSNV_RF_SCORE	Interpro_domain	AAChange.knownGene	rmsk	MetaSVM_score	Freq_gnomAD_genome_POPs	OMIM	Phenotype_MIM	OrphaNumber	Orpha	Otherinfo
# 1	925952	925952	G	A	SAMD11	exonic	nonsynonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.G11A:p.G4E,SAMD11:ENST00000616016.5:exon2:c.G548A:p.G183E,SAMD11:ENST00000618323.5:exon2:c.G548A:p.G183E	SAMD11:NM_152486:exon2:c.G11A:p.G4E	clinvar: Uncertain_significance 	 InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[1, 0, 0, 0, 0, 0, 0, 0] 	.	.	.	4.311	29.7	0.0	3.6	.	.	.	SAND domain|SAND domain|SAND domain;.;.;.;.;.;.;.;.;.;.	SAMD11:ENST00000616125.5:exon1:c.G11A:p.G4E,SAMD11:ENST00000617307.5:exon1:c.G11A:p.G4E,SAMD11:ENST00000618181.5:exon1:c.G11A:p.G4E,SAMD11:ENST00000618779.5:exon1:c.G11A:p.G4E,SAMD11:ENST00000622503.5:exon1:c.G11A:p.G4E,SAMD11:ENST00000342066.8:exon2:c.G11A:p.G4E,SAMD11:ENST00000437963.5:exon2:c.G11A:p.G4E,SAMD11:ENST00000616016.5:exon2:c.G548A:p.G183E,SAMD11:ENST00000618323.5:exon2:c.G548A:p.G183E	.	-0.858	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	925956	925956	C	T	SAMD11	exonic	synonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.C15T:p.I5I,SAMD11:ENST00000616016.5:exon2:c.C552T:p.I184I,SAMD11:ENST00000618323.5:exon2:c.C552T:p.I184I	SAMD11:NM_152486:exon2:c.C15T:p.I5I	clinvar: UNK 	 InterVar: Likely benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 1, 0, 0, 1, 0] 	.	.	.	.	.	.	.	.	.	.	.	SAMD11:ENST00000616125.5:exon1:c.C15T:p.I5I,SAMD11:ENST00000617307.5:exon1:c.C15T:p.I5I,SAMD11:ENST00000618181.5:exon1:c.C15T:p.I5I,SAMD11:ENST00000618779.5:exon1:c.C15T:p.I5I,SAMD11:ENST00000622503.5:exon1:c.C15T:p.I5I,SAMD11:ENST00000342066.8:exon2:c.C15T:p.I5I,SAMD11:ENST00000437963.5:exon2:c.C15T:p.I5I,SAMD11:ENST00000616016.5:exon2:c.C552T:p.I184I,SAMD11:ENST00000618323.5:exon2:c.C552T:p.I184I	.	.	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	925969	925969	C	T	SAMD11	exonic	nonsynonymous SNV	SAMD11	rs200686669	SAMD11:ENST00000342066.8:exon2:c.C28T:p.P10S,SAMD11:ENST00000616016.5:exon2:c.C565T:p.P189S,SAMD11:ENST00000618323.5:exon2:c.C565T:p.P189S	SAMD11:NM_152486:exon2:c.C28T:p.P10S	clinvar: UNK 	 InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[1, 0, 0, 0, 0, 0, 0, 0] 	0.0004	.	.	4.121	27.9	0.0	4.52	.	.	.	.;.;.;.;.;.;.;.;.;.;.	SAMD11:ENST00000616125.5:exon1:c.C28T:p.P10S,SAMD11:ENST00000617307.5:exon1:c.C28T:p.P10S,SAMD11:ENST00000618181.5:exon1:c.C28T:p.P10S,SAMD11:ENST00000618779.5:exon1:c.C28T:p.P10S,SAMD11:ENST00000622503.5:exon1:c.C28T:p.P10S,SAMD11:ENST00000342066.8:exon2:c.C28T:p.P10S,SAMD11:ENST00000437963.5:exon2:c.C28T:p.P10S,SAMD11:ENST00000616016.5:exon2:c.C565T:p.P189S,SAMD11:ENST00000618323.5:exon2:c.C565T:p.P189S	.	-0.469	AFR:0.0001,AMR:0,EAS:0,FIN:0.0020,NFE:0.0002,OTH:0.0020,ASJ:0	616765	.			.
# 1	925976	925976	T	C	SAMD11	exonic	nonsynonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.T35C:p.I12T,SAMD11:ENST00000616016.5:exon2:c.T572C:p.I191T,SAMD11:ENST00000618323.5:exon2:c.T572C:p.I191T	SAMD11:NM_152486:exon2:c.T35C:p.I12T	clinvar: UNK 	 InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[1, 0, 0, 0, 0, 0, 0, 0] 	.	.	.	4.213	28.8	0.0	4.52	.	.	.	.;.;.;.;.;.;.;.;.;.;.	SAMD11:ENST00000616125.5:exon1:c.T35C:p.I12T,SAMD11:ENST00000617307.5:exon1:c.T35C:p.I12T,SAMD11:ENST00000618181.5:exon1:c.T35C:p.I12T,SAMD11:ENST00000618779.5:exon1:c.T35C:p.I12T,SAMD11:ENST00000622503.5:exon1:c.T35C:p.I12T,SAMD11:ENST00000342066.8:exon2:c.T35C:p.I12T,SAMD11:ENST00000437963.5:exon2:c.T35C:p.I12T,SAMD11:ENST00000616016.5:exon2:c.T572C:p.I191T,SAMD11:ENST00000618323.5:exon2:c.T572C:p.I191T	.	-0.553	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	925986	925986	C	T	SAMD11	exonic	synonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.C45T:p.C15C,SAMD11:ENST00000616016.5:exon2:c.C582T:p.C194C,SAMD11:ENST00000618323.5:exon2:c.C582T:p.C194C	SAMD11:NM_152486:exon2:c.C45T:p.C15C	clinvar: UNK 	 InterVar: Likely benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 1, 0, 0, 1, 0] 	.	.	.	.	.	.	.	.	.	.	.	SAMD11:ENST00000616125.5:exon1:c.C45T:p.C15C,SAMD11:ENST00000617307.5:exon1:c.C45T:p.C15C,SAMD11:ENST00000618181.5:exon1:c.C45T:p.C15C,SAMD11:ENST00000618779.5:exon1:c.C45T:p.C15C,SAMD11:ENST00000622503.5:exon1:c.C45T:p.C15C,SAMD11:ENST00000342066.8:exon2:c.C45T:p.C15C,SAMD11:ENST00000437963.5:exon2:c.C45T:p.C15C,SAMD11:ENST00000616016.5:exon2:c.C582T:p.C194C,SAMD11:ENST00000618323.5:exon2:c.C582T:p.C194C	.	.	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	926003	926003	C	T	SAMD11	exonic	nonsynonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.C62T:p.S21F,SAMD11:ENST00000616016.5:exon2:c.C599T:p.S200F,SAMD11:ENST00000618323.5:exon2:c.C599T:p.S200F	SAMD11:NM_152486:exon2:c.C62T:p.S21F	clinvar: UNK 	 InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[1, 0, 0, 0, 0, 0, 0, 0] 	.	.	.	4.149	28.2	0.001	4.63	.	.	.	.;.;.;.;.;.;.;.;.;.;.	SAMD11:ENST00000616125.5:exon1:c.C62T:p.S21F,SAMD11:ENST00000617307.5:exon1:c.C62T:p.S21F,SAMD11:ENST00000618181.5:exon1:c.C62T:p.S21F,SAMD11:ENST00000618779.5:exon1:c.C62T:p.S21F,SAMD11:ENST00000622503.5:exon1:c.C62T:p.S21F,SAMD11:ENST00000342066.8:exon2:c.C62T:p.S21F,SAMD11:ENST00000437963.5:exon2:c.C62T:p.S21F,SAMD11:ENST00000616016.5:exon2:c.C599T:p.S200F,SAMD11:ENST00000618323.5:exon2:c.C599T:p.S200F	.	-0.468	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	926010	926010	G	T	SAMD11	exonic	synonymous SNV	SAMD11	.	SAMD11:ENST00000342066.8:exon2:c.G69T:p.P23P,SAMD11:ENST00000616016.5:exon2:c.G606T:p.P202P,SAMD11:ENST00000618323.5:exon2:c.G606T:p.P202P	SAMD11:NM_152486:exon2:c.G69T:p.P23P	clinvar: UNK 	 InterVar: Likely benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 1, 0, 0, 1, 0] 	.	.	.	.	.	.	.	.	.	.	.	SAMD11:ENST00000616125.5:exon1:c.G69T:p.P23P,SAMD11:ENST00000617307.5:exon1:c.G69T:p.P23P,SAMD11:ENST00000618181.5:exon1:c.G69T:p.P23P,SAMD11:ENST00000618779.5:exon1:c.G69T:p.P23P,SAMD11:ENST00000622503.5:exon1:c.G69T:p.P23P,SAMD11:ENST00000342066.8:exon2:c.G69T:p.P23P,SAMD11:ENST00000437963.5:exon2:c.G69T:p.P23P,SAMD11:ENST00000616016.5:exon2:c.G606T:p.P202P,SAMD11:ENST00000618323.5:exon2:c.G606T:p.P202P	.	.	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.
# 1	926014	926014	G	A	SAMD11	splicing	.	SAMD11	.	.	.	clinvar: UNK 	 InterVar: Pathogenic PVS1=1 PS=[0, 0, 0, 0, 0] PM=[0, 1, 0, 0, 0, 0, 0] PP=[0, 0, 1, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0] 	.	.	.	5.758	34	.	4.63	.	1.0000	0.938	.;.;.;.;.;.;.;.;.;.;.	.	.	.	AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.	616765	.			.

vcf_head = """##fileformat=VCFv4.2
##INFO=<ID=INTERVAR_CLASS,Number=.,Type=String,Description="InterVar overall classification">
##INFO=<ID=INTERVAR_PVS1,Number=.,Type=Integer,Description="InterVar PVS1 score">
##INFO=<ID=INTERVAR_PS,Number=.,Type=String,Description="InterVar PS score matrix">
##INFO=<ID=INTERVAR_PM,Number=.,Type=String,Description="InterVar PM score matrix">
##INFO=<ID=INTERVAR_PP,Number=.,Type=String,Description="InterVar PP score matrix">
##INFO=<ID=INTERVAR_BA1,Number=.,Type=String,Description="InterVar BA1 score matrix">
##INFO=<ID=INTERVAR_BS,Number=.,Type=String,Description="InterVar BS score matrix">
##INFO=<ID=INTERVAR_BP,Number=.,Type=String,Description="InterVar BP score matrix">
#CHROM\tPOS\tID	REF\tALT\tQUAL\tFILTER\tINFO
"""
in_file = args.intervar_tsv
out_file = args.output_basename + ".vcf.gz"
# open bgzip pipe for compression
p = subprocess.Popen("bgzip -c > " + out_file, shell=True, stdin=subprocess.PIPE)
p.stdin.write(bytes(vcf_head, 'utf-8'))
tsv = open(in_file)
skip = next(tsv)
for line in tsv:
    info = line.rstrip('\n').split("\t")
    record = "chr" + info[0] + "\t" + info[1] + "\t.\t" + info[3] + "\t" + info[4]
    intervar_info = info[13]
    # parse this absolute gem of a string:
    #  InterVar: Pathogenic PVS1=1 PS=[0, 0, 0, 0, 0] PM=[0, 1, 0, 0, 0, 0, 0] PP=[0, 0, 1, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0] 
    parsed = re.match(" InterVar: (.*)\s(PVS1=\d)\s(PS=.*)\s(PM=.*)\s(PP=.*)\s(BA1=.*)\s(BS=.*)\s(BP=.*)\s", intervar_info)
    record += "\t.\t.\tINTERVAR_CLASS=" +  ";".join(parsed.groups()) +  "\n"
    p.stdin.write(bytes(record, 'utf-8'))
p.communicate()

subprocess.call("tabix " + out_file, shell=True)