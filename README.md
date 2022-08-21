# MSGENOVA: Variants detection and analyse_mutation_site

[![Snakemake](https://img.shields.io/badge/snakemake-≥ 7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

bla bla bla

## Authors

* Gabrièle Adam (gabriele.adam@inrae.fr)
* Nadia Bessoltane (nadia.bessoltane@inrae.fr)
* Delphine Charif (delphine.charif@inrae.fr)

## 1. How to install the pipeline on Genotoul
This pipeline has been develloped to run on the Genotoul cluster (http://bioinfo.genotoul.fr), and the workflow management is done with snakemake. The user must have a Genotoul account (http://bioinfo.genotoul.fr/index.php/ask-for/create-an-account-2/). 
### 1.1. For all users (excluding IJPB users)
The pipeline has already been installed in the spae reserved for the IJPB institute. To run the pipeline, you need to be given access to the code. Please contact the bioinformatic team at bioinfo-ijpb@inrae.fr.
### 1.2. For IJPB users
For IJPB members, the pipeline is already install on the cluster in /save/projects/ijpb/bioinfo-code/src/msgenova folder

## 2. Usage
Once the pipeline is installed, fill the configuration file, prepare the correct input file, and finally, run the pipeline.
### 2.1. Fill the configuration file with pathways for the differents input
The configuration file contains pathways to all the inputs needed to run the pipeline. It's a json formatted file, meaning, a line begins by an attribute, next a colon (":"), and next a value for this attribute. It must be created by the user in order to run the pipeline. The config file for the pipeline must contain at least 4 attributes:
* **sample***: pathway to the csv formateed table containing  barcodes informations
* **reference**: pathway to the reference, that must be in fasta format
* **gff**: pathway to the annotation file in gff format
* **project_name**: name for the project
* **DP.min: Minimum sequencing depth for variant filtering
* **AR.min: Min AR value for variant filtering
```bash
#Exemple of configuration file
sample: "pathway_to_sample_sheet"
reference: "pathway_to_reference"
gff: "pathway_to_reference_annotation"
project_name: "project_name"
AR.min: 0.4
DP: 2
```
sample, reference, gff and project name are mandatory attributes. Other attributes can also be specified:
* **regions**: 4 columns file containing coordinates to the regions of interests (optional). This file should have 4 columns: chromosome, begining of the region, end of the region, name of the region
* **vector** : Vector sequences in fasta format (in order to look for vector insertion in the data) (optional)
* **min.DP** : Minimum sequencing depth for variant analysis
* **min.AR**: Minimum AD value for variant analysis
* **cpu**: cpu given for each job of the analysis (default: 8)
* **ram**: ram given for each job of the pipeline (default: 30G)
### 2.2 Input files
### 2.2.1. Sample sheet (mandatory)
The sample sheet is a tab separated 4 columns files. The 4 columns are named :
* **sample**: sample name, 
* **condition**: condition for this sample. There must always be at least one reference sample named wt in the sheet, 
* **fq1**: Pathway to the forward read for this sample, 
* **fq2**: Pathway to the reverse read for this sample
```bash
#example of sample sheet
sample	condition	fq1	fq2
S1	wt	S1.R1.fq.gz	S1.R2.fq.gz
S2	mt	S2.R1.fq.gz	S2.R2.fq.gz
```
The condition column must always contains at least one value **wt**. This reference will be used for variant filtering.
### 2.2.2. Reference (mandatory)
Fasta formatted genome of reference
### 2.2.3. Annotation file (mandatory)
Annotation for the reference in gff3 format
```bash
#example of gff2 file for arabidopsis thaliana
1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;Name=NAC001;biotype=protein_coding;description=NAC domain-containing protein 1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q0WV96];gene_id=AT1G01010;logic_name=araport11
```
### 2.2.4. Project name (mandatory)
Name for the project. It is better for it not to have any tabulation or whitespace.
### 2.2.5. Regions file (optional)
3 columns tabulation separated file containing the regions specific on which the analysis will be done. The first columns should be the chromosome name, the second columns should be the beginning of the region, and the third column should be the end of the region.
```bash
1	10	100
3	500	1000
```
### 2.2.6. Vector file (optional)
Fasta formatted file containing the different vector sequences. This file will be used to look for plasmidic insertion in the genome using the tool tDNAscan (https://doi.org/10.3389/fgene.2019.00685). The script was modified to be implemented in our pipeline.
### 2.2.7. Minimal DP (optional)
Minimal sequencing depth value used to filter the variant during the analysis
### 2.2.8. Minimal AD (optional)
Minimal number of sample carying the mutation. Value used to filter variant during quality control.
### 2.2.9. cpu (optional)
Number of cpu used to run each job (default: 8)
### 2.2.10. ram (optional)
Size of the ram assigned to each job (default: 30G)
### 2.3 Quickstart
To run the pipeline, once the configuration file is filled, please use the command:
```bash
./submit -c pathway_to_the_configuration_file -o pathway_to_the_output_directory
```
The 2 arguments are the configuration file (filled as mentionned above), and the output directory. This directory must exists **before** the pipeline is run.

## 3. Output files
The output files are present in the results folder that will be created by the pipeline. The different folder in the results directory are:
```bash
├── 00_logs
├── 01_sequence_qc
├── 02_demultiplex
├── 03_mapping
├── 04_snv_indels
├── 05_sv_calling
├── 06_haplotype_pbaa
├── genome
``` 
**00_logs**  contains the log file for each job send on the cluster
**01_sequence_qc** contains the multiqc results in a html report called Project_name_multiqc_trim.html.
**02_mapping** contains 4 subfolder:
 - bam: results of the alignment
 - coverage: results of the coverage done with samtools coverage
 - depth: results of the depth analysis done with samtools
 - flagstat: results of the flagstat command with flagstat
**03_snv_indels_calling**: contains results of the variant calling done with haplotypecaller. Files are named Project_name.snpeff.html (with results from snpeff analysis), and Project_name_snv_indel.vcf (vcf containing the variant found during the analysis)
**04**: 
**05_tdnascan**: contains result of the vector search on the data. The search is done by sample eand by vector present in the vector file.
**genome**: contains result of the reference indexation

## 4. Versions of the tools used in the pipeline
Several tools are used in this workflow:
- bwa
- gatk4 version 4.2.6.1
- picard 
- bedtools 
- samtools
- fastqc
- multiqc
- trimmomatic
- R version 4.0.4
- SnpEff version 5.0
- pandoc version 2.1.3
- ucsc-fasplit

## 5. FAQ
### Q1. Ma job was killed on the cluster / I killed my job
In case your job was killed/ you killed your job, the pipeline need to be rerun. Beforehand, due to snakemake way of working, the folder .snakemake/locks in the working directory need to be deleted with the command rm.
```bash
# Example of command line ton delete the folder
rm -R .snakemake/locks
```
### Q2. I want to install this pipeline on another cluster than Genotoul
This pipeline can be installed on another cluster than Genotoul. To be able to run the pipeline on another system, the git repository need to by clone first, and second, the configuration file need to be filled.
#### 1. Set up your working environement
To be able to run the pipeline, the user must have access to the forgemia (https://forgemia.inra.fr/), and be a member of the project, as the repos is still private. Next, to be able to clone the repos, the user must create a couple of ssh keys following the instruction present on forgemia (https://forgemia.inra.fr/-/profile/keys). The workflow can then, be installed on the cluster.
#### 2. Install snakemake 
This pipeline has been coded using snakemake version 7.8.1. This is the minimal version needed to run the pipeline. If snakemake is not install in your environment, you can sinply install it with conda/mamba by using these commands:
```bash
#create the snakemake environment
mamba create -c conda-forge -c bioconda -n snakemake snakemake
#activate the environment
conda activate snakemake
```
You can also contact your cluster administration team to install it for you.
#### 3. How to install the pipeline on a computer cluster
Several options are available to download the pipeline from git. 
##### 3.1 Download the tarball
You can either download the source code in zip, tar.gz, tar.bz2 or tar format. In that case, please decompress the folder in your working directory
##### 3.2 Clone the repos
Else, you can clone the git repository with ssh or https protocol with these commands:
```bash
#command to clone git repository using ssh
git clone git@forgemia.inra.fr:sps-bioinfo/multicrispr.git [folder_name(optionnal)]
#or using https
git clone https://forgemia.inra.fr/sps-bioinfo/multicrispr.git [folder_name(optionnal)]
```
The pipeline is now downloaded on the cluster. 
##### 4.1 Fill the configuration file
As indicated above
##### 4.2 Fill the advanced configuration file
The first configuration file must contains the pathway to the files needed to run the pipeline. The second configuration file contains the options for several tools, and also the pathway to tools used by the pipeline but installed locally on the Genotoul cluster. PLease change these pathways, for the local installation on your own cluster system.
##### 5. Run  the pipeline
As explained above
### Q3. What to do if the pipeline stopped before the end ?
If a problem arised and the pipeline stop, please delete the .snakemake/locks folder in your working directory and rerun the pipeline
### Q4. How is the pipeline run on a cluster
This script launch the snakemake workflow on z cluster using a sbatch command. Options -p -c --mem are the ones from the cluster (p= partition on which the job is sent, c=number of CPUs to use to run the job, mem=RAM allocated to the job). Other options are the ones from snakemake:
* **snakefile**: Path to Snakefile file containing the rule "all"
<!--  * **use-conda**: Use conda to load and create environment according to yaml files in envs/n -->
* **p**: Print command line
* **config-file**: Rule configuration
* **cluster-config**: Path to cluster configuration for each job
* **cluster**: Use cluster mode of snakemake 
* **jobs**: Number of jobs launched at the same time by snakemake
* **latency-wait**: Time that snakemake will wait in seconds after the job finishes, allowing snakemake to check if output files are present
* **use-conda**: Use conda environment to run each job 
* **stats**: print statistics on the time taken by each job
```bash
#Example of submit_job file 
#!/bin/bash
config="pathway_to_configuration_file"
path="pathway_to_working_directory"
sbatch snakemake -p --snakefile $path/Snakefile --configfile $config --cluster-config $path/cluster_config.json --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.mem} --output {cluster.output} --error {cluster.error} --job-name {rule} --mem {resources.mem_mb} --cpu-per-tasks {threads}" --jobs 100 --latency-wait 90 --use-envmodules --stats stat.json --directory $path  >& snakemake.log
```
Once you have a stable release for your snakemake workflow, you may export a rule graph showing the logical execution flow chart.  
  ```bash
snakemake --rulegraph | dot -Tpng >your-snakemake-workflow_rulegraph.png
```  
Graphviz library allows you to have different output format ex:pdf, png,jpg. To obtain another format, please modify the -T option and the extension of the output rulegraph file.
### Q5. What is the workflow management used to develop this pipeline
This pipeline was coded using snakemake version 7.8.1
### Q6. How is organized the git repository ?
The basic workflow project tree is as follows:  
```bash
# inside the cloned snakemake workflow template directory
├── cluster_config.json
├── config.yaml
├── config_advanced.yaml
├── README.md
├── rules/
├── script/
├── Snakefile
├── submit_job
```
* `config.yaml`: snakemake workflow configuration file  containing pathway to the different files needed to run the pipeline
* `config_advanced.yaml`: advanced snakemake workflow configuration file containing options for several job and pathway to tool used by the pipeline
* `LICENSE`: license file (should be open source: CC BY-SA 4.0)   
* `README.md`: this file (to be renamed or deleted)    
* `Snakefile`: main snakemake workflow file   
* `cluster_config.json`: cluster configuration in json format for each rule and by default
* `submit_job`: submission file
* `rules`: snakemake workflow subrules files (\*.smk)   
* `script`: local script neede to create the final report (\*.Rmd)   
### Q7. How to define the different parameters of a job send on the cluster ?
The parameters of each job send on the cluster are defined in the file cluster_config.json. Advanced user can modify these parameters if needed.
```bash
#Example of a default job configuration
    "__default__" :
    {   
        "partition" : "All",
        "c" : "16",
        "mem":1G,
        "name" : "{rule}",
        "output" : "logs/cluster/{rule}.{wildcards.plant}.out",
        "error" : "logs/cluster/{rule}.{wildcards.plant}.err",
        "N" : "1"
    }
```
You can fill the configuration information for all rules. In the example above, the slurm parameters used are:
- partition: partition on which to send the job
- c : number of CPUs on which to run the job
- mem: RAM allocated to the job
- name: name given to the job
- output: output file of the job
- error: error file of the job
- N: number of cores on which to run the job


# LICENSE  
  
CC BY-SA 4.0  
    
