rule create_input_pindel:
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
		wd=WORKDIR
	conda:
		"../envs/pindel.yaml"
	shell:
		"""
		pindel -T {threads} -f {input.reference} -i {input.config_pindel} -o {params.wd}/results/04_pindel/{wildcards.sample}
		"""
