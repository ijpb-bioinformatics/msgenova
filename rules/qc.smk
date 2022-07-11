rule multiqc:
	input:
		expand("results/01_sequence_qc/trimmed.{s.sample}.{R}_fastqc.{ext}",ext=["zip"],R=["R1","R2"],s=SAMPLE.itertuples()),
	output:
		"results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
#	log:
#		err=WORKDIR+"/essai_1320.err",
#		out=WORKDIR+"/essai_1320.out"
	params:
		extra=config["params_multiqc"]
	conda:
		"../envs/qc.yaml"
	shell:
		"""
		directory=`dirname {input[0]}`
		directory_res=`dirname {output}`
		multiqc --force -n {output} {input}
		"""
#	wrapper:
#		"v1.7.0/bio/multiqc"

rule fqc_essai2:
	input:
		"/save/project/ijpb/bioinfo-data/raw-data/SOJ/fastq/JAP3-jmj14-4_TAATGCGC-AGGCTATA_L007_R2.fastq.gz"
	output:
		html="qc/JAP3-jmj14-4_TAATGCGC-AGGCTATA_L007_R2_fastqc.html",
		zip="qc/JAP3-jmj14-4_TAATGCGC-AGGCTATA_L007_R2_fastqc.zip"
	params: " --quiet"
	threads: get_thread
	message: "Do fqc2"
	log:
		"log/essai_log_file.log"
	resources:
		mem_mb=config["ram"]
	shell:
		"""
		echo "$(date '+%Y%m%d %r') [$(basename $0): QC & Trimming] Starting Sequence Quality Control and Trimming Process..." | tee -a {log} 2>&1
		fastqc {params} -t {threads} --outdir /work/gadam/Mudetect/qc {input}  
		"""

rule trimming:
	input:
		r1=get_fq1,
		r2=get_fq2
	output:
		r1="results/01_sequence_qc/trimmed.{sample}.R1.fq.gz",
		r2="results/01_sequence_qc/trimmed.{sample}.R2.fq.gz",
		r1_unpaired="results/01_sequence_qc/trimmed.{sample}.R1.unpaired.fq.gz",
		r2_unpaired="results/01_sequence_qc/trimmed.{sample}.R2.unpaired.fq.gz",
		log="results/01_sequence_qc/{sample}.trimmomatic.log"
	params:
		trimmer=["TRAILING:3"],
		extra=config["params_trimmomatic"],
		compression_level="-9"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	conda:
		"../envs/qc.yaml"
#	log:
#		error="results/01_sequence_qc/log/{sample}.trimmomatic.error",
#		output="results/01_sequence_qc/log/{sample}.trimmomatic.output"
	shell:
		"""
		trimmomatic PE -threads {threads} {params.extra} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} {params.trimmer} 2> {output.log}
		"""


rule trimmed_fqc2:
	input:
		"results/01_sequence_qc/trimmed.{sample}.R2.fq.gz",
	output:
		html="results/01_sequence_qc/trimmed.{sample}.R2_fastqc.html",
		zip="results/01_sequence_qc/trimmed.{sample}.R2_fastqc.zip"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	conda:
		"../envs/qc.yaml"
	params:
		extra=config["params_fastqc"]
	shell:
		"""
		fastqc --outdir results/01_sequence_qc/ {params.extra} {input}
		"""

rule trimmed_fqc1:
	input:
		"results/01_sequence_qc/trimmed.{sample}.R1.fq.gz",
	output:
		html="results/01_sequence_qc/trimmed.{sample}.R1_fastqc.html",
		zip="results/01_sequence_qc/trimmed.{sample}.R1_fastqc.zip"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	conda:
		"../envs/qc.yaml"
	params:
		extra=config["params_fastqc"]
	shell:
		"""
		fastqc --outdir	results/01_sequence_qc/	{params.extra} {input}
		"""

