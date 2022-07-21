"""
Author: Gabrièle Adam, Nadia Bessoltane, Delphine Cherif
Affiliation: IJPB
Aim : Analysis of variants from Illumina data

Date: 14th June 2022
"""
import re
from re import search
import pandas as pd
import sys
import glob
from snakemake.utils import validate, min_version
import itertools
import subprocess
from pathlib import Path
from os.path import exists

##### set minimum snakemake version #####
min_version("5.32.0")

## Define workdir 
WORKDIR=os.getcwd()
configfile: config["repo_script"]+"/config_advanced"
include: "rules/common.smk"

SAMPLE = pd.read_table(config["sample"], dtype=str,delimiter="\t").set_index(["sample"], drop=False)
#SAMPLE.index = SAMPLE.index.set_levels([i.astype(str) for i in SAMPLE.index.levels]) 

#create results and logs output directory
LIST_of_DIR=["results/00_logs/","results/01_sequence_qc/log/","results/03_snv_indels_calling/log/","results/02_mapping/log/","results/06_report/log/"]
for dir in LIST_of_DIR:
	if not (os.path.isdir(dir)):
		create_logsdir(dir)

# Define reference name 
REF_NAME=Path(config["reference"]).stem

#Define reference variable (in case of concatenation of reference with vector)
REFERENCE=define_reference(config,"vector")

print(REFERENCE)

#Define bwa suffix 
SUFFIX_BWA=["amb","ann","bwt","pac","sa"]
QC_LIST=["depth","coverage","flagstat"]

LIST_REGION_HAPLOTYPECALLER=[]
def read_reference(reference):
	LIST_CHR_REFERENCE=[]
	dict={}
	i=0
	# try to get the chromosome from fasta
	with open(reference,"r") as file:
	#with open(config["reference"],"r") as file:
		for line in file:
			if re.search(">",line):
				#if ">" is present then cut by space, and get first sequence
				split_line=line.split(" ")
				result=re.sub(">","",split_line[0])
				#print(result)
				LIST_CHR_REFERENCE.append(i)
				dict[str(i)]=re.sub(r"[\n\t\s]*", "", result)
				i=i+1
				#LIST_CHR_REFERENCE.append(re.sub(r"[\n\t\s]*", "", result))
	return LIST_CHR_REFERENCE,dict


#regions: "/shared/projects/camelina_pacbio_pipe/msgenova/regions.bed"

#Define bool_region for booleen value if there is a region file present
bool_region = False

#We are creating a dictionary containing number as keys and chromosome or pathway to region file as value
dict_chr={}

# essai fct for optional config
#def get_config_key(cfg,key,default_list=None):
def get_config_key(cfg,key,bool_var,dict):
	try:
		value=cfg[key]
		print("regions file exist")
		bool_var= True
		#Define Regions for analysis
		#REGIONS = pd.read_table(config["regions"],dtype=str,delimiter="\t",header=None,names=["chromosome","beg","end"]).set_index(["chromosome","beg"], drop=False)
		LIST_CHR=[]
		dict={}
		dict[str(0)]=config["regions"]
		LIST_CHR.append(0)
		#for elm in REGIONS.itertuples():
		#	LIST_CHR.append(elm.chromosome+":"+elm.beg+"-"+elm.end)
		return LIST_CHR,dict
		#return value
	except:
		print("no regions file - Read reference")
		bool_var= False
		default_list=read_reference(config["reference"])
		return default_list

LIST_REGION_HAPLOTYPECALLER,dict_chr=get_config_key(config,"regions",bool_region,dict_chr)
print("List regions haplotypecaller")
print(LIST_REGION_HAPLOTYPECALLER,dict_chr.values(),dict_chr.keys())

#essai dictionnaire
x=dict_chr.get(0)
print(x)
x=dict_chr.get(1)
print(x)
#Transform project name
NAME_PROJECT=transform_project_name(config["project_name"])

rule all:
	input:
		#"results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html",
		##expand("results/02_mapping/bam/{s.sample}.bam",s=SAMPLE.itertuples()),
		#expand("results/01_sequence_qc/{s.sample}.trimmomatic.log",s=SAMPLE.itertuples()),
		#expand("results/02_mapping/{qc}/{s.sample}.{qc}", qc=QC_LIST,s=SAMPLE.itertuples()),
		#"results/03_snv_indels_calling/variants.merged.vcf",
		#"results/03_snv_indels_calling/"+NAME_PROJECT+".html"
		#"results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html",
		#"results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf"
		#"results/06_report/msgenova_report.html"
		"results/genome/"+REFERENCE+".fasta"

include: "rules/qc.smk"
include: "rules/index.smk"
include: "rules/align.smk"
include: "rules/qc_bam.smk"
include: "rules/sv_calling.smk"
include: "rules/annotate_variant.smk"
include: "rules/report.smk"
