rule build_snpeff_db:
	"""
	To be able to annotate the variant using snpeff, a database need to be created prior. It contains all variant already present in the genome of reference.
	"""
	output:
		config="results/genome/snpeff.config"
	params:
		ref_name=REFERENCE,
		ver_snpeff="SnpEff-5.0",
		ref="results/genome/"+REFERENCE+".fasta",
		gff=config["gff"],
		wk=WORKDIR,
		extra=config["params_snpeff_build"],
		vector=get_vector(config,"vector")[0],
		ref_config=config["reference"]
	conda:
		"../envs/snpeff.yaml"
	resources:
		mem_mb=get_mem
	threads: get_thread	
	shell:
		"""
		echo {params.ref_name}.genome > {output.config}
		if [ ! -d results/genome/{params.ref_name} ]
		then
			mkdir -p results/genome/{params.ref_name}
		fi
		echo -e "{params.ref_name}\n{params.ver_snpeff}" > {params.ref_name}.db
		if [ ! -L {params.wk}/results/genome/{params.ref_name}/sequences.fa ]
		then
			if [ {params.vector} == "TRUE" ]
			then
				ln -s {params.wk}/{params.ref} {params.wk}/results/genome/{params.ref_name}/sequences.fa
			else
				ln -s {params.ref_config} {params.wk}/results/genome/{params.ref_name}/sequences.fa
			fi
		fi
		if [ ! -L {params.wk}/results/genome/{params.ref_name}/genes.gff ]
		then
			ln -s {params.gff} {params.wk}/results/genome/{params.ref_name}/genes.gff
		fi
		snpEff build -c {output.config} {params.extra} -gff3 -v {params.ref_name} -dataDir .
		"""

rule annotate_variant:
	"""
	Snpeff is used to annotate variant. The output is a html file containing analysis of the variants obtained.
	"""
	input:
		vcf="results/03_snv_indels_calling/variants.merged.vcf",
		config="results/genome/snpeff.config"
	output:
		vcf="results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf",
		html="results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html"
	params:
		ref_name=REFERENCE,
		extra=config["params_snpeff_ann"]
	conda:
		"../envs/snpeff.yaml"
	resources:
		mem_mb=get_mem
	threads: get_thread
	shell:
		"""
		snpEff eff -c {input.config} {params.extra} -dataDir . {params.ref_name} -s {output.html} {input.vcf} > {output.vcf}
		"""

