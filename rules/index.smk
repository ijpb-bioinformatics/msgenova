rule bwa_index_reference:
	"""
	Index reference for bwa
	"""
	input:
		ref=config["reference"]
	output:
		index=expand("results/genome/"+REF_NAME+".fasta.{suffix}",suffix=SUFFIX_BWA),
		ref="results/genome/"+REF_NAME+".fasta"
	conda:
		"../envs/align.yaml"
	resources:
		mem_mb=config["ram"]
	threads: get_thread
	params:
		extra=config["params_bwa_index_reference"]
	shell:
		"""
		ln -s {input.ref} {output.ref}
		bwa index {params.extra} {output.ref}
		"""

rule create_dict_reference:
	"""
	To be able to run snpeff, a dictionnary need to be created using picard
	"""
	input:
		ref="results/genome/"+REF_NAME+".fasta"
	output:
		ref="results/genome/"+REF_NAME+".dict"
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_create_dict_reference"]
	shell:
		"""
		picard -Xmx{resources.mem_mb} CreateSequenceDictionary --REFERENCE {input.ref} --OUTPUT {output.ref} {params.extra}
		"""

rule samtools_index_reference:
	input:
		ref="results/genome/"+REF_NAME+".fasta"
	output:
		ref="results/genome/"+REF_NAME+".fasta.fai"
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_index_reference_samtools"]
	shell:
		"""
		samtools faidx {params.extra} {input.ref}
		"""


rule index_alignment_file:
	"""
	Index alignment file so that HaploypeCaller could be run
	"""
	input:
		bam="results/02_mapping/bam/{sample}.bam"
	output:
		bai="results/02_mapping/bam/{sample}.bam.bai"
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	params:
		extra=config["params_index_alignment"]
	shell:
		"""
		samtools index -@ {threads} {params.extra} {input.bam}
		"""

#rule index_intersect_alignment_file:
#	input:
#		bam="results/tmp/03_mapping/by_file/{well}.intersect.bam"
#	output:
#		bai="results/tmp/03_mapping/by_file/{well}.intersect.bam.bai"
#	conda:
#		"../envs/picard.yaml"
#	threads: get_thread
#	resources:
#		mem_mb=config["ram"]
#	shell:
#		"""
#		samtools index -@ {threads} {input.bam}
#		"""
#
