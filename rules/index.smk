rule deal_with_vector:
	"""
	Concatenate reference with vector if it is present in the config, else, just symlink to reference
	"""
	input:
		ref=config["reference"]
	output:
		ref="results/genome/"+REFERENCE+".fasta"
	params:
		vector=get_vector(config,"vector")[0],
		field=get_vector(config,"vector")[1]
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"""
		val_bool={params.vector}
		if [ $val_bool == "TRUE" ]
		then
			cat {input.ref} {params.field} > {output.ref}
		else
			ln -s {input.ref} {output.ref}
		fi
		"""

rule bwa_index_reference:
	"""
	Index reference for bwa
	"""
	input:
		ref="results/genome/"+REFERENCE+".fasta"
		#ref=get_reference
		#ref=config["reference"]
	output:
		index=expand("results/genome/"+REFERENCE+".fasta.{suffix}",suffix=SUFFIX_BWA),
		#ref="results/genome/"+REFERENCE+".fasta"
	conda:
		"../envs/tdnascan.yaml"
	resources:
		mem_mb=get_mem
		#mem_mb=config["ram"]
	threads: get_thread
	params:
		extra=config["params_bwa_index_reference"]
	shell:
		"""
		bwa index {params.extra} {input.ref}
		"""

rule create_dict_reference:
	"""
	To be able to run snpeff, a dictionnary need to be created using picard
	"""
	input:
		ref="results/genome/"+REFERENCE+".fasta"
	output:
		ref="results/genome/"+REFERENCE+".dict"
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		extra=config["params_create_dict_reference"]
	shell:
		"""
		picard -Xmx{resources.mem_mb} CreateSequenceDictionary --REFERENCE {input.ref} --OUTPUT {output.ref} {params.extra}
		"""

rule samtools_index_reference:
	input:
		ref="results/genome/"+REFERENCE+".fasta"
	output:
		ref="results/genome/"+REFERENCE+".fasta.fai"
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
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
		mem_mb=get_mem
	params:
		extra=config["params_index_alignment"]
	shell:
		"""
		samtools index -@ {threads} {params.extra} {input.bam}
		"""
