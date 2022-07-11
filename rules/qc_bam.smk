rule flagstat:
	input:
		"results/02_mapping/bam/{sample}.mark_duplicates.bam"
	output:
		"results/02_mapping/flagstat/{sample}.flagstat"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_flagstat"]
	shell:
		"""
		samtools flagstat {params.extra} {input} > {output}
		"""

rule samtools_coverage:
	input:
		"results/02_mapping/bam/{sample}.bam"
	output:
		"results/02_mapping/coverage/{sample}.coverage"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_samtools_coverage"]
	shell:
		"""
		samtools coverage -o {output} {params.extra} {input}
		"""

rule samtools_depth:
	input:
		"results/02_mapping/bam/{sample}.bam"
	output:
		"results/02_mapping/depth/{sample}.depth"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_samtools_depth"]
	shell:
		"""
		samtools depth {params.extra} {input} > {output}
		"""
