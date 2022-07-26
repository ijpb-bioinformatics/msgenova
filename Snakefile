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
LIST_of_DIR=["results/00_logs/","results/01_sequence_qc/log/","results/02_mapping/log/"]
for dir in LIST_of_DIR:
	if not (os.path.isdir(dir)):
		create_logsdir(dir)

# Define reference name
REF_NAME=Path(config["reference"]).stem

#Define reference variable (in case of concatenation of reference with vector)
REFERENCE=define_reference(config,"vector")

#vector: "/shared/projects/camelina_pacbio_pipe/msgenove_21/vector.fa"

#if vector
if get_vector(config,"regions")[0] == "TRUE":
	print("there is a region file")
	print(get_vector(config,"regions")[0])
	REGIONS=pd.read_table(get_vector(config,"regions")[1], dtype=str,delimiter="\t",header=None,names=["chr","beg","end"]).set_index(["chr","beg"], drop=False)

#Define bwa suffix
SUFFIX_BWA=["amb","ann","bwt","pac","sa"]
QC_LIST=["depth","coverage","flagstat"]

LIST_REGION_HAPLOTYPECALLER=[]
#regions: "/shared/projects/camelina_pacbio_pipe/msgenova/regions.bed"

#Define bool_region for booleen value if there is a region file present
bool_region = False

#We are creating a dictionary containing number as keys and chromosome or pathway to region file as value
dict_chr={}

LIST_REGION_HAPLOTYPECALLER,dict_chr=get_config_key(config,"regions",bool_region,dict_chr)
print("List regions haplotypecaller")
print(LIST_REGION_HAPLOTYPECALLER,dict_chr.values(),dict_chr.keys())

#Transform project name
NAME_PROJECT=transform_project_name(config["project_name"])

#define list of vector name if vector is present in the configuration file
#if get_vector(config,"vector")[0] == "TRUE":
#	LIST_VECTOR,dict_vector=read_reference(get_vector(config,"vector")[1])

def get_all(wildcards):
	if get_vector(config,"vector")[0] == "TRUE":
		LIST_VECTOR,dict_vector=read_reference(get_vector(config,"vector")[1])
		expand("results/05_tdnascan/{s.sample}/{vector}/5.{vector}_insertion.bed",s=SAMPLE.itertuples(),vector=LIST_VECTOR),
rule all:
	input:
		#expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),
		#expand("results/02_mapping/coverage/{s.sample}_{reg.chr}:{reg.beg}-{reg.end}.coverage",s=SAMPLE.itertuples(),reg=REGIONS.itertuples()),
		#expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		#"results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html",
		#"results/06_report/msgenova_report.html"
		#expand("results/05_tdnascan/{s.sample}/5.tdna_insertion.bed",s=SAMPLE.itertuples()),
		#get_all
		"aggregate.txt"

include: "rules/qc.smk"
include: "rules/index.smk"
include: "rules/align.smk"
include: "rules/qc_bam.smk"
include: "rules/sv_calling.smk"
include: "rules/annotate_variant.smk"
include: "rules/report.smk"
include: "rules/tdnascan.smk"
