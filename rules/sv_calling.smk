def get_region(wildcards):
	x= dict_chr.get(str(wildcards.regions))
	#print(x)
	return x

rule HaplotyCaller:
	"""
	HaplotypeCaller from the GATK package calls SNP and indels by sample, using a local de-novo assembly of haplotypes. When a variation is observed in a region, the reads mapped to this region are re-assembled in that region. The mode used here is GVCF. Samples will be next genotype together using GenotypeGVCFs from GATK. 
	"""
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		ref="results/genome/"+REFERENCE+".fasta",
		index="results/genome/"+REFERENCE+".fasta.fai",
		dict="results/genome/"+REFERENCE+".dict",
		bai="results/02_mapping/bam/{sample}.bam.bai"
	output:
		vcf=temp("results/03_snv_indels_calling/{sample}.{regions}.g.vcf"),
		idx=temp("results/03_snv_indels_calling/{sample}.{regions}.g.vcf.idx")
	conda:
		"../envs/env_gatk4.yaml"
	threads: get_thread("HaplotyCaller")
	params:
		extra=config["params_haplotype_caller"],
		#regions= lambda wildcards,output: get_region(wildcards)
		#region=dict_chr.get(wildcards.regions)
		regions=get_region
	resources:
		mem_mb=get_mem("HaplotyCaller")
	shell:
		"""
		 gatk --java-options -Xmx{resources.mem_mb} HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -L {params.regions} --native-pair-hmm-threads {threads} {params.extra}
		"""

rule GenomicsDB:
	"""
	Here GenomicsDB replace CombineGVCFs because it was taking too much time to run. GenomicsDB create a database and  store all variations present in gvcf file. Later sample will be genotyped
	"""
	input:
		vcf=expand("results/03_snv_indels_calling/{s.sample}.{{regions}}.g.vcf",s=SAMPLE.itertuples()),
		idx=expand("results/03_snv_indels_calling/{s.sample}.{{regions}}.g.vcf.idx",s=SAMPLE.itertuples()),
	output:
		dir=temp(directory("results/03_snv_indels_calling/genomicDB_{regions}")),
		json="results/03_snv_indels_calling/genomicDB_{regions}/callset.json",
		vcf="results/03_snv_indels_calling/genomicDB_{regions}/vcfheader.vcf",
		vidmap="results/03_snv_indels_calling/genomicDB_{regions}/vidmap.json"
	params:
		region=get_region,
		extra=config["params_genomics_db"]
	conda:
		"../envs/env_gatk4.yaml"
	threads: get_thread("GenomicsDB")
	resources: mem_mb=get_mem("GenomicsDB")
	shell:
		"""
		if [ -d {output.dir} ]
		then
			rm -R {output.dir}
		fi
		variant=""
		for file in {input.vcf}
		do
			variant=$variant" -V "$file
		done
		gatk --java-options -Xmx{resources.mem_mb} GenomicsDBImport $variant --genomicsdb-workspace-path {output.dir} --intervals {params.region} {params.extra}
		"""

rule Genotype_gvcf:
	"""
	Genotype samples using gatk GenotypeGVCFs
	"""
	input:
		ref="results/genome/"+REFERENCE+".fasta",
		dir="results/03_snv_indels_calling/genomicDB_{regions}/",
		#db=expand("results/03_snv_indels_calling/genomicDB/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"]),
		db=expand("results/03_snv_indels_calling/genomicDB_{{regions}}/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"]),
	output:
		vcf=temp("results/03_snv_indels_calling/variants.{regions}.vcf"),
		idx=temp("results/03_snv_indels_calling/variants.{regions}.vcf.idx"),
	conda:
		"../envs/env_gatk4.yaml"
	resources: mem_mb=get_mem("Genotype_gvcf")
	threads: get_thread("Genotype_gvcf")
	params:
		extra=config["params_genotype_gvcf"]
	shell:
		"""
		list_db=""
		for file in {input.dir}
		do
			list_db=$list_db" gendb://"$file" "
		done
		echo $list_db
		gatk GenotypeGVCFs --reference {input.ref} -V gendb://{input.dir} -O {output.vcf} {params.extra}
		"""

rule gatherVCF:
	input:
		vcf=expand("results/03_snv_indels_calling/variants.{regions}.vcf",regions=LIST_REGION_HAPLOTYPECALLER),
		idx=expand("results/03_snv_indels_calling/variants.{regions}.vcf.idx",regions=LIST_REGION_HAPLOTYPECALLER),
	output:
		vcf=temp("results/03_snv_indels_calling/variants.merged.vcf"),
		idx=temp("results/03_snv_indels_calling/variants.merged.vcf.idx")
	conda:
		"../envs/env_alignment.yaml"
	resources: mem_mb=get_mem("gatherVCF")
	threads: get_thread("gatherVCF")
	params:
		extra=config["params_gathervcf"]
	shell:
		"""
		list=""
		for file in {input.vcf}
		do
			list=$list" I="$file
		done
		picard -Xmx{resources.mem_mb} GatherVcfs $list O={output.vcf} {params.extra}
		"""
