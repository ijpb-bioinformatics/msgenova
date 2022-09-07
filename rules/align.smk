rule align:
	"""
	Align data using bwa mem
	"""
	input:
		fq1="results/01_sequence_qc/trimmed.{sample}.R1.fq.gz",
		fq2="results/01_sequence_qc/trimmed.{sample}.R2.fq.gz",
		index=expand("results/genome/"+REFERENCE+".fasta.{suffix}",suffix=SUFFIX_BWA),
		reference="results/genome/"+REFERENCE+".fasta"
	output:
		sam=temp("results/02_mapping/{sample}.sam")
	params:
		extra=config["params_bwa"]
	conda:
		"../envs/tdnascan.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"""
		bwa mem -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:Illumina\\tPL:Illumina\\tPU:{wildcards.sample}\" -M -t {threads} {params.extra} {input.reference} {input.fq1} {input.fq2} > {output.sam} 
		"""

rule sort_sam:
	"""
	Use picard to sort sam in alignment
	"""
	input:
		"results/02_mapping/{sample}.sam"
	output:
		bam=temp("results/02_mapping/{sample}.sort.bam"),
		bai=temp("results/02_mapping/{sample}.sort.bai")
	params:
		sort_order="coordinate",
		extra=config["params_sort_sam"]
	threads: get_thread
	resources:
		mem_mb=get_mem
	conda:
		"../envs/picard.yaml"
	shell:
		"""
		picard SortSam  -Xmx{resources.mem_mb} --INPUT {input} --OUTPUT {output.bam} --SORT_ORDER {params.sort_order} {params.extra}
		"""

rule mark_duplicate:
	"""
	Use picard to mark duplicates
	"""
	input:
		"results/02_mapping/{sample}.sort.bam"
	output:
		#enlever temp pour avoir les vecteurs
		sam=temp("results/02_mapping/{sample}.mark_duplicates.bam"),
		metrics="results/02_mapping/{sample}.mark_duplicates.metrics"
	threads: get_thread
	resources:
		mem_mb=get_mem
	conda:
		"../envs/picard.yaml"
	params:
		extra=config["params_mark_duplicate"]
	shell:
		"""
		picard MarkDuplicates -Xmx{resources.mem_mb} --INPUT {input} --OUTPUT {output.sam} -M {output.metrics} {params.extra}
		"""


#rule RealignerTargetCreator:
#	"""
#	Realign Indels using gatk
#	"""
#	input:
#		reference="results/genome/"+REF_NAME+".fasta",
#		sam="results/02_mapping/{sample}.{condition}.sort.sam"
#	output:
#		"results/02_mapping/{sample}.{condition}.realign_indels.sam"
#	threads: get_thread
#	resources:
#		mem_mb=config["ram"]
#	conda:
#		"../envs/gatk3.yaml"
#	shell:
#		"""
#		java -Xmx{resources.mem_mb} -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R {input.reference} -I {input.sam} -o {output} -rf BadCigar --num_threads {threads} -fixMisencodedQuals
#		""" 


rule keep_mapped_only:
	"""
	Keep mapped reads only in alignment and intersect if a region file is present
	"""
	input:
		"results/02_mapping/{sample}.mark_duplicates.bam"
	output:
		"results/02_mapping/bam/{sample}.bam"
	resources:
		mem_mb=get_mem
	conda:
		"../envs/picard.yaml"
	threads: get_thread
	params:
		extra=config["params_extra_mapped_reads"],
		bool=get_vector(config,"regions")[0],
		region=get_vector(config,"regions")[1]
	shell:
		"""
		
		if [ {params.bool} == "TRUE" ]
		then
			samtools view --threads {threads} {params.extra} -F4 -h -S -b {input} | bedtools intersect -a - -b {params.region} > {output}
		else
			samtools view --threads {threads} {params.extra} -F4 -h -S -b {input} > {output}
		fi
		"""
