def get_path_vector(wildcards):
	return get_vector(config,"vector")[1]

checkpoint cut_vector_file:
	input:
		vector=get_path_vector
	output:
		directory("results/genome/tdnascan/vectors/")
	conda:
		"../envs/tdnascan.yaml"
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

rule copy_tdna:
	input:
		vector=get_path_vector
	output:
		tdna_seq="results/genome/TDNA_sequence.fasta"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wk=WORKDIR
	shell:
		"""
		cp {input} {params.wk}/{output.tdna_seq}
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
		"../envs/tdnascan.yaml"
	params:
		extra=config["params_bwa_index_reference"],
		wk=WORKDIR,
	shell:
		"""
		bwa index {params.extra} {input.vector}
		"""

rule prepare_reference_tdnascan:
	"""
	Prepare symbolic link for reference
	"""
	input:
		ref=config["reference"]
		#ref=get_vector(config,"reference")[1]
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
		"../envs/tdnascan.yaml"
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
		bed="results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.bed",
		sam=temp("results/05_tdnascan/{sample}/{vector}/1.TDNA.sam"),
		bam1=temp("results/05_tdnascan/{sample}/{vector}/1.TDNA_sort.bam"),
		sam_insert=temp("results/05_tdnascan/{sample}/{vector}/5.insertion.sam")
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
		path_script=`readlink -f {params.install_dir}/script/tdnascan.py`
		cd results/05_tdnascan/{wildcards.sample}
		# add argument for install parameters
		python2.7 $path_script -1 {params.wk}/{input.fq1} -2 {params.wk}/{input.fq2} -t {params.wk}/{input.vector} -g {params.wk}/{input.reference} -p {wildcards.vector} -@ {threads} -i {params.wk}/script -d {params.wk}/results/05_tdnascan/{wildcards.sample} {params.extra}
		"""

rule clean_and_delete_tdnascan:
	input:
		sam="results/05_tdnascan/{sample}/{vector}/1.TDNA.sam",
		bam1="results/05_tdnascan/{sample}/{vector}/1.TDNA_sort.bam",
		sam_insert="results/05_tdnascan/{sample}/{vector}/5.insertion.sam",
		bed="results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.bed"
	output:
		"results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.reduce.bed"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wd=WORKDIR
	shell:
		"""
			rm -R {params.wd}/results/05_tdnascan/{wildcards.sample}/{wildcards.vector}/tmp
			awk ' {{if ( $6 > 0 ) print $0}}' {input.bed} > {output}
		"""

rule tdnascan_annotate:
	input:
		"results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.reduce.bed"
	output:
		"results/05_tdnascan/{sample}/{vector}/5.{vector}_insertion.annotated.new.bed"
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		gff=get_vector(config,"gff")[1],
		install_dir=config["repo_script"],
		wd=WORKDIR
	conda:
		"../envs/tdnascan.yaml"
	shell:
		"""
		path_script=`readlink -f {params.install_dir}/script/tdnaAnnot.py`
		cd results/05_tdnascan/{wildcards.sample}
		python2.7 $path_script -i {params.wd}/{input} -f {params.gff} -o {params.wd}/{output}
		"""



#def aggregate_vector(wildcards):
#	#list=[]
#	checkpoint_output=checkpoints.cut_vector_file.get(**wildcards).output[0]
 #	return expand("results/05_tdnascan/{s.sample}/{v}/5.{v}_insertion.annotated.new.bed", s=SAMPLE.itertuples(),v=glob_wildcards(os.path.join(checkpoint_output, "{v}.fa")).v)
#	#list.extend(expand("results/05_tdnascan/{s.sample}/{v}/5.{v}_insertion.annotated.new.bed", s=SAMPLE.itertuples(),v=glob_wildcards(os.path.join(checkpoint_output, "{v}.fa")).v))
#	#return list
#

#rule aggregate:
#	input:
#		aggregate_vector
#	output:
#		"aggregate.txt"
#	shell:
#		"""
#		for file in {input}
#		do
#			echo $file >> {output}
#		done
#		"""
#
#
