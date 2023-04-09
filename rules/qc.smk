rule multiqc:
	input:
		expand("results/01_sequence_qc/trimmed.{s.sample}.{R}_fastqc.{ext}",ext=["zip"],R=["R1","R2"],s=SAMPLE.itertuples()),
	output:
		"results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html"
	threads: get_thread("multiqc")
	resources:
		mem_mb=get_mem("multiqc")
	params:
		extra=config["params_multiqc"]
	conda:
		"../envs/env_qc.yaml"
	shell:
		"""
		directory=`dirname {input[0]}`
		directory_res=`dirname {output}`
		multiqc --force -n {output} {input}
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
		log="results/01_sequence_qc/log/{sample}.trimmomatic.log"
	params:
		trimmer=["TRAILING:3"],
		extra=config["params_trimmomatic"],
		compression_level="-9"
	threads: get_thread("trimming")
	resources:
		mem_mb=get_mem("trimming")
	conda:
		"../envs/env_qc.yaml"
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
	threads: get_thread("trimmed_fqc2")
	resources:
		mem_mb=get_mem("trimmed_fqc2")
	conda:
		"../envs/env_qc.yaml"
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
	threads: get_thread("trimmed_fqc1")
	resources: mem_mb=get_mem("trimmed_fqc1")
	conda:
		"../envs/env_qc.yaml"
	params:
		extra=config["params_fastqc"]
	shell:
		"""
		fastqc --outdir	results/01_sequence_qc/	{params.extra} {input}
		"""

rule concatenate_log_trimmomatic:
	input:
		expand("results/01_sequence_qc/log/{s.sample}.trimmomatic.log", s=SAMPLE.itertuples()),
	output:
		"results/01_sequence_qc/log/trimmomatic.log"
	threads: get_thread("concatenate_log_trimmomatic")
	resources: mem_mb=get_mem("concatenate_log_trimmomatic")
	params:
		wd=WORKDIR
	shell:
		"""
		line=`cat {input[0]} | grep Surviving| sed 's/)/)\\n/g' | sed 's/Both/\\nBoth/g' | sed 's/:.*//g' | sed 's/$/-/g'| sed 's/\t//g' | sed 's/ //g'| tr -d "\n" | sed 's/-/\t/g' `
		echo -e "Sample\t"$line >> {output}
		for file in {input}
		do
			s_name=`basename {params.wd}/$file | sed 's/\.trimmomatic\.log//g' `
			line=`cat $file | grep Surviving| sed 's/)/)\\n/g' | sed 's/Both/\\nBoth/g' | sed 's/.*://g' | sed 's/$/-/g'| sed 's/\t//g' | sed 's/ //g'| tr -d "\n" | sed 's/-/\t/g' `	
			echo -e $s_name"\t"$line >> {output}
		done
		"""
