# MSGENOVA: Multi-sample comprehensive analysis of genomic variations in Illumina DNA libraries. 

[![Snakemake](https://img.shields.io/badge/snakemake-≥ 7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

bla bla bla

## Authors

* Gabrièle Adam (gabriele.adam@inrae.fr)
* Nadia Bessoltane (nadia.bessoltane@inrae.fr)
* Delphine Charif (delphine.charif@inrae.fr)

## 1. Installation 
The pipeline has benn tested on Genotoul and IFB computing clusters. It can run on a slurm computing system, and on a local computer (in that case, the coputing time will be long)
For IJPB members, the pipeline is already installed in the ijpb project space on the genotoul server. Please contact the bioinformatic team at bioinfo-ijpb@inrae.fr.
### 1.2. Prerequisite
This pipeline has been develloped to run on a slurm cluster, with the workflow management tool **snakemake (minimal version 7.7.0)** and **conda (minimum version 4.4.10)**.
#### 1.2.1. Install conda and snakemake
To install conda, please follow the instructions on https://conda.io/projects/conda/en/latest/user-guide/install/index.html
If snakemake is not install in your environment, you can sinply install it with conda/mamba and activate the environment by using these commands:
```bash
#create the snakemake environment
mamba create -c conda-forge -c bioconda -n snakemake snakemake
#activate the environment
conda activate snakemake
```
You can also contact your cluster administration team to install it for you.
### 1.3. Download the source code
The user must have access to the forgemia (https://forgemia.inra.fr/), and be a member of the project, as the repos is still private. Next, to be able to clone the repos, the user must create a couple of ssh keys following the instruction present on forgemia (https://forgemia.inra.fr/-/profile/keys). The workflow can then, be installed on the cluster. The code can be downloaded or clone via git.
##### 1.2.2.1 Download the tarball
Download the tarball using the wget command:
```bash
wget https://forgemia.inra.fr/sps-bioinfo/msgenova/-/archive/main/msgenova-main.tar.gz
```
You can download the source code in either zip, tar.gz, tar.bz2 or tar format. Next, please decompress the folder in your working directory
##### 1.2.2.2 Clone the repos
You can clone the git repository with ssh or https protocol with these commands (the folder_name option allows the user to use the name of the new directory):
```bash
#command to clone git repository using ssh
git clone git@forgemia.inra.fr:sps-bioinfo/msgenova.git [folder_name(optionnal)]
#or using https
git clone https://forgemia.inra.fr/sps-bioinfo/msgenova.git [folder_name(optionnal)]
```
The pipeline is now downloaded on the cluster. 
### 1.4. Configuration
Options configurations are present in the config_advanced file. They must be only modify by an advanced user.
### 1.5. Test

```bash

###

```

## 2. Usage
Once the pipeline is installed, fill the configuration file, prepare the correct input file, and finally, run the pipeline.

To run the pipeline, once the configuration file is filled, please use the command:
```bash
./msgenova -c pathway_to_the_configuration_file.txt -o pathway_to_the_output_directory [-s] [-r]
```

| Option | Description |
|--------|----------------------------------------------------------------------------------------------------------|
| -c | configuration file (filled as indicated above) |
| -o | output directory (must exist before running the pipeline) |
| -s | run indels search, run by default if a vector file is present |
| -r | AR value for variant filtering FALSE by default |

### User configuration file 
The configuration file contains pathways to all the inputs needed to run the pipeline. It's a json formatted file, meaning, a line begins by an attribute, next a colon (":"), and next a value for this attribute. It must be created by the user. Arguments are listed below, and the mandatory one are in bold.


| attribute        | Description  |
|------------------|------------------------------------------------------------------------------------|
| **project_name** | name for the project. It is better for it not to have any tabulation or whitespace |
| **sample**       | pathway to the sample sheet(see a)        |
| **reference**    | pathway to the reference, that must be in fasta format  |
| **gff**          | pathway to the annotation file in gff format  |
| DP.min           | Minimum sequencing depth for variant filtering  |
| AR.min           | Min AR value for variant filtering  |
| regions          | 4 columns file containing coordinates to the regions of interests (optional). This file should have 4 columns: chromosome, begining of the region, end of the region, name of the region (see b)  |
| vector           | Vector sequences in fasta format (in order to look for vector insertion in the data) (see c.) |
| cpu              | cpu given for each job of the analysis (default: 8)  |
| mem              | ram given for each job of the pipeline (default: 30G)  |


```bash
#Exemple of configuration file
project_name: "project_name"
sample: "pathway_to_sample_sheet"
reference: "pathway_to_reference"
gff: "pathway_to_reference_annotation"
AR.min: 0.4
DP: 2
```

#### a. Sample sheet

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

#### b. Regions file 
3 columns tabulation separated file containing the regions specific on which the analysis will be done. The first columns should be the chromosome name, the second columns should be the beginning of the region, and the third column should be the end of the region.
```bash
1	10	100
3	500	1000
```
#### c. Vector file 
Fasta formatted file containing the different vector sequences. This file will be used to look for plasmidic insertion in the genome using the tool tDNAscan (https://doi.org/10.3389/fgene.2019.00685). The script was modified to be implemented in our pipeline.

## 3. Output files
The output files are present in the **results** folder that will be created by the pipeline. The different sub folder present are:
```bash
├── 00_logs
├── 01_sequence_qc
├── 02_mapping
├── 03_snv_indels_calling
├── 04_pindel
├── 05_tdnascan
├── 06_report
├── genome
``` 
* **00_logs**  contains the log file for each job send on the cluster
* **01_sequence_qc** contains the multiqc results in a html report called Project_name_multiqc_trim.html.
* **02_mapping** contains 4 subfolders:
* - bam: results of the alignment
* - coverage: results of the coverage done with samtools coverage
* - depth: results of the depth analysis done with samtools
* - flagstat: results of the flagstat command with flagstat
* **03_snv_indels_calling**: contains results of the variant calling done with haplotypecaller. Files are named Project_name.snpeff.html (with results from snpeff analysis), and Project_name_snv_indel.vcf (vcf containing the variant found during the analysis)
* **04_pindel**: contains results from pindel tool. Results are present in the file Project_nqme_sv.vcf.gz.
* **05_tdnascan**: contains result of the vector search on the data. The search is done by sample eand by vector present in the vector file.
* **genome**: contains reference and vectors index

## 4. Versions of the tools used in the pipeline
Several tools are used in this workflow:
- bwa version 0.7.12
- gatk4 version 4.2.6.1
- picard version 2.27.4
- bedtools version 2.30.0
- samtools version 1.15.1
- tabix version 1.11
- bcftools version 1.15.1
- python version 3 and 2.7
- fastqc version 0.11.9
- multiqc version 1.13a
- trimmomatic version 0.39
- SnpEff version 5.0
- pandoc version 2.1.3
- ucsc-fasplit version 377
- pindel version 0.2.5b9
- R version 4.1.3

## 5. FAQ
### Q1. Ma job was killed on the cluster / I killed my job
In case your job was killed/ you killed your job, the pipeline need to be rerun. Beforehand, due to snakemake way of working, the folder .snakemake/locks in the working directory need to be deleted with the command rm.
```bash
# Example of command line ton delete the folder
rm -R .snakemake/locks
```
### Q2. What to do if the pipeline stopped before the end ?
If a problem arised and the pipeline stop, please delete the .snakemake/locks folder in your working directory and rerun the pipeline
### Q3. How is the pipeline run on a cluster
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
### Q4. How is organized the git repository ?
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
### Q5. How to define the different parameters of a job send on the cluster ?
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
    
