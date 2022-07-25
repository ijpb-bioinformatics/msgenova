rule flagstat:
	input:
		"results/02_mapping/{sample}.mark_duplicates.bam"
	output:
		"results/02_mapping/flagstat/{sample}.flagstat"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		extra=config["params_flagstat"]
	shell:
		"""
		samtools flagstat {params.extra} {input} > {output}
		"""

rule samtools_coverage_by_regions:
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		bai="results/02_mapping/bam/{sample}.bam.bai"
	output:
		"results/02_mapping/coverage/{sample}_{chr}:{beg}-{end}.coverage"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		extra=config["params_samtools_coverage"]
	shell:
		"""
			samtools coverage -r {wildcards.chr}:{wildcards.beg}-{wildcards.end} -o {output} {params.extra} {input.bam}
		"""

def get_input_samtools_coverage(wildcards):
	if get_vector(config,"regions")[0] == "TRUE":
		LIST=[]
		REGIONS=pd.read_table(get_vector(config,"regions")[1], dtype=str,delimiter="\t",header=None,names=["chr","beg","end"]).set_index(["chr","beg"], drop=False)
		for elm in REGIONS.itertuples():
			#print(elm)
			LIST.append("results/02_mapping/coverage/"+wildcards.sample+"_"+elm.chr+":"+elm.beg+"-"+elm.end+".coverage")
		return LIST
	else:
		return "results/02_mapping/coverage/{sample}.temp.coverage"

rule samtools_coverage:
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		bai="results/02_mapping/bam/{sample}.bam.bai"
	output:
		"results/02_mapping/coverage/{sample}.temp.coverage"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		extra=config["params_samtools_coverage"],
		bool=get_vector(config,"regions")[0],
		regions=get_vector(config,"regions")[1]
	shell:
		"""
		samtools coverage -o {output} {params.extra} {input.bam}
		"""

rule samtools_coverage_final:
	input:
		get_input_samtools_coverage
	output:
		"results/02_mapping/coverage/{sample}.coverage"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		bool=get_vector(config,"regions")[0],
	shell:
		"""
		if [ {params.bool} == "TRUE" ]
		then
			cat {input[0]} | grep "#" > {output}
			cat {input} | grep -v "#" >> {output}
		else
			cp {input} {output}
		fi
		"""
		
#samtools coverage -r {params.regions} -o {output} {params.extra} {input.bam}
rule samtools_depth:
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		bai="results/02_mapping/bam/{sample}.bam.bai"
	output:
		"results/02_mapping/depth/{sample}.depth"
	conda:
		"../envs/samtools.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		extra=config["params_samtools_depth"],
		bool=get_vector(config,"regions")[0],
		regions=get_vector(config,"regions")[1]
	shell:
		"""
		if [ {params.bool} == "TRUE" ]
		then
			samtools depth -b {params.regions} {params.extra} {input.bam} > {output}
		else
			samtools depth {params.extra} {input.bam} > {output}
		fi
		"""
