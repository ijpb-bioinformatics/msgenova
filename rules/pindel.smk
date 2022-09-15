rule create_input_pindel:
	"""
	Create configuration file for pindel
	"""
	input:
		"results/02_mapping/bam/{sample}.bam"
	output:
		"results/04_pindel/input_pindel_{sample}.txt"
	params:
		wd=WORKDIR,
		ins_size=get_size_insert
	threads: get_thread
	resources: mem_mb=get_mem
	shell:
		"""
		if [ ! -d {params.wd}/results/04_pindel/ ]
		then
			mkdir -p {params.wd}/results/04_pindel/
		fi
		nom=`basename {input}`
		echo -e "{params.wd}/{input}\t{params.ins_size}\t$nom" >> {output}
		"""

rule run_pindel:
	"""
	Run pindel sample by sample
	"""
	input:
		config_pindel="results/04_pindel/input_pindel_{sample}.txt",
		files="results/02_mapping/bam/{sample}.bam",
		reference="results/genome/"+REFERENCE+".fasta",
		index_bam="results/02_mapping/bam/{sample}.bam.bai",
		fai_ref="results/genome/"+REFERENCE+".fasta.fai"
	output:
		expand("results/04_pindel/{{sample}}_{ext_pindel}",ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wd=WORKDIR,
		extra=config["params_run_pindel"]
	conda:
		"../envs/pindel.yaml"
	shell:
		"""
		pindel {params.extra} -T {threads} -f {input.reference} -i {input.config_pindel} -o {params.wd}/results/04_pindel/{wildcards.sample}
		"""

rule convert_pindel:
	"""
	Convert pindel output to vcf format
	"""
	input:
		expand("results/04_pindel/{{sample}}_{ext_pindel}",ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),
	output:
		temp("results/04_pindel/{sample}.vcf")
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		reference="results/genome/"+REFERENCE+".fasta",
		wd=WORKDIR,
		extra=config["params_convert_pindel"]
	conda:
		"../envs/pindel.yaml"
	shell:
		"""
		pindel2vcf {params.extra} -r {params.reference} -P {params.wd}/results/04_pindel/{wildcards.sample} -R xxx -d xx -v {params.wd}/results/04_pindel/{wildcards.sample}.vcf
		"""


rule merge_vcf_files:
	"""
	Merge all samples together
	"""
	input:
		expand("results/04_pindel/{s.sample}.vcf", s=SAMPLE.itertuples()),
	output:
		temp("results/04_pindel/pindel_results.vcf.gz")
	threads: get_thread
	resources: mem_mb=get_mem
	conda:
		"../envs/picard.yaml"
	params:
		extra=config["params_merge_vcf_files"],
		wd=WORKDIR
	shell:
		"""
		list=""
		i=0
		for file in {input}
		do
			if [[ -f {params.wd}/$file.gz ]]
			then
				rm {params.wd}/$file.gz
			fi
			if [[ -f {params.wd}/$file.gz.tbi ]]
			then
				rm {params.wd}/$file.gz.tbi
			fi
			bgzip {params.wd}/$file
			tabix {params.wd}/$file.gz
			list=$list" "{params.wd}"/"$file".gz"
			((i=i+1))
		done
		#if there are several files in in the input then merge, else only move because bcftools merge do not work if there is only one file
		if [ $i == "1" ]
		then
			cp $list {output}
		else
			bcftools merge -m none -Oz {params.extra} -o {output} $list
		fi
		"""

rule annotate_pindel_output:
	"""
	Annotate pindel putput (in vcf format) using Snpeff
	"""
	input:
		vcf="results/04_pindel/pindel_results.vcf.gz",
		config="results/genome/snpeff.config"
	output:
		vcf="results/04_pindel/"+NAME_PROJECT+"_sv.vcf.gz",
		html="results/04_pindel/"+NAME_PROJECT+"_sv.html"
	params:
		ref_name=REFERENCE,
		extra=config["params_snpeff_ann"]
	conda:
		"../envs/snpeff.yaml"
	resources: mem_mb=get_mem
	threads: get_thread
	shell:
		"""
		snpEff -Xmx{resources.mem_mb} eff -c {input.config} {params.extra} -dataDir . {params.ref_name} -s {output.html} {input.vcf} > {output.vcf}
		"""
