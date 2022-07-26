checkpoint cut_vector_file:
	input:
		vector=get_vector(config,"vector")[1]
	output:
		directory("results/genome/tdnascan/vectors/")
		#expand("results/genome/tdnascan/vectors/{vector}.fa",vector=LIST_VECTOR)
	conda:
		"../envs/ucsc.yaml"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wk=WORKDIR
	shell:
		"""
		mkdir -p {params.wk}/results/genome/tdnascan/
		mkdir -p {params.wk}/results/genome/tdnascan/vectors/
		faSplit byname {input} {params.wk}/results/genome/tdnascan/vectors/
		"""

rule index_vector:
	"""
	Index vector for tdnascan
	"""
	input:
		vector="results/genome/tdnascan/vectors/{vector}.fa"
	output:
		index=expand("results/genome/tdnascan/vectors/{{vector}}.fa.{suffix}",suffix=SUFFIX_BWA),
	threads: get_thread
	resources: mem_mb=get_mem
	conda:
		"../envs/align.yaml"
	params:
		extra=config["params_bwa_index_reference"],
		wk=WORKDIR,
	shell:
		"""
		bwa index {params.extra} {input.vector}
		"""


#		cd {params.wk}
#		if [ -d {params.wk}/05_tdnascan/{wildcards.sample}/{wildcards.vector} ]
#		then
#			rm -R {params.wk}/05_tdnascan/{wildcards.sample}/{wildcards.vector}
#		fi
#		#create output directory in results/05_tdnascan
#		mkdir -p {params.wk}/results/05_tdnascan/{wildcards.sample}/{wildcards.vector}
#		#In this dir put a symlink for vector
#		ln -s {params.wk}/{input.vector} {params.wk}/results/05_tdnascan/{wildcards.sample}/{wildcards.vector}/{wildcards.vector}.fa
#		#now index this vector in the results/genome directory
#		bwa index {params.extra} {params.wk}/results/05_tdnascan/{wildcards.sample}/{wildcards.vector}/{wildcards.vector}.fa
#		# create  symlink to 
#		"""
#

rule prepare_reference_tdnascan:
	"""
	Prepare symbolic link for reference
	"""
	input:
		ref=get_vector(config,"reference")[1]
	output:
		ref="results/genome/tdnascan/reference/"+REF_NAME+".fa",
		index=expand("results/genome/tdnascan/reference/"+REF_NAME+".fa.{suffix}",suffix=SUFFIX_BWA),	
		#"results/05_tdnascan/{wildcards.sample}/{wildcards.vector}/"+REF_NAME+".fa"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wk=WORKDIR,
		extra=config["params_bwa_index_reference"]
	conda:
		"../envs/align.yaml"
	shell:
		"""
		ln -s {input.ref} {params.wk}/{output.ref}
		bwa index {params.extra} {params.wk}/{output.ref}
		"""

rule tdnascan:
	input:
		fq1="results/01_sequence_qc/trimmed.{sample}.R1.fq.gz",
		fq2="results/01_sequence_qc/trimmed.{sample}.R2.fq.gz",
		reference="results/genome/tdnascan/reference/"+REF_NAME+".fa",
		vector="results/genome/tdnascan/vectors/{vector}.fa",
		index_ref=expand("results/genome/tdnascan/reference/"+REF_NAME+".fa.{suffix}",suffix=SUFFIX_BWA),
		index_vector=expand("results/genome/tdnascan/vectors/{{vector}}.fa.{suffix}",suffix=SUFFIX_BWA),
	output:
		"results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.bed"
	params:
		install_dir=config["repo_script"],
		wk=WORKDIR,
		extra=config["params_tdnascan"]
	threads: get_thread
	resources: mem_mb=get_mem
	conda:
		"../envs/tdnascan.yaml"
	shell:
		"""
		cd {params.wk}
		if [ -d {params.wk}/05_tdnascan/{wildcards.sample}/{wildcards.vector} ]
		then
			rm -R {params.wk}/05_tdnascan/{wildcards.sample}/{wildcards.vector}
		fi
		mkdir -p results/05_tdnascan/{wildcards.sample}/{wildcards.vector}
		cd results/05_tdnascan/{wildcards.sample}
		# add argument for install parameters
		python2.7 {params.install_dir}/script/tdnascan.py -1 {params.wk}/{input.fq1} -2 {params.wk}/{input.fq2} -t {params.wk}/{input.vector} -g {params.wk}/{input.reference} -p {wildcards.vector} -@ {threads} -i {params.wk}/script -d {params.wk}/results/05_tdnascan/{wildcards.sample} {params.extra}
		"""

def aggregate_vector(wildcards):
	checkpoint_output=checkpoints.cut_vector_file.get(**wildcards).output[0]
	return expand("results/05_tdnascan/{s.sample}/{v}/5.{v}_insertion.bed", s=SAMPLE.itertuples(),v=glob_wildcards(os.path.join(checkpoint_output, "{v}.fa")).v)
 
rule aggregate:
	input:
		aggregate_vector
	output:
		"aggregate.txt"
	shell:
		"""
		for file in {input}
		do
			echo $file >> {output}
		done
		"""
