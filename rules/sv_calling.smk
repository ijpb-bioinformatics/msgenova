def get_region(wildcards):
	x= dict_chr.get(str(wildcards.regions))
	#print("essai dictionaire")
	#print(wildcards.regions)
	#for key in dict_chr.keys():
	#	print("Key :"+str(key))
	#for val in dict_chr.values():
	#	print("Value :"+str(val))
	#print("test keys")
	#if wildcards.regions in dict_chr.keys():
	#	print("in the list")
	#else:
	#	print("not in the list")
	print(x)
	return x

rule HaplotyCaller:
	"""
	HaplotypeCaller from the GATK package calls SNP and indels by sample, using a local de-novo assembly of haplotypes. When a variation is observed in a region, the reads mapped to this region are re-assembled in that region. The mode used here is GVCF. Samples will be next genotype together using GenotypeGVCFs from GATK. 
	"""
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		ref="results/genome/"+REF_NAME+".fasta",
		index="results/genome/"+REF_NAME+".fasta.fai",
		dict="results/genome/"+REF_NAME+".dict",
		bai="results/02_mapping/bam/{sample}.bam.bai"
	output:
		vcf=temp("results/03_snv_indels_calling/{sample}.{regions}.g.vcf")
	conda:
		"../envs/gatk4.yaml"
	threads: get_thread
	params:
		extra=config["params_haplotype_caller"],
		#regions= lambda wildcards,output: get_region(wildcards)
		#region=dict_chr.get(wildcards.regions)
		regions=get_region
	resources:
		mem_mb=config["ram"]
	shell:
		"""
		 gatk --java-options -Xmx{resources.mem_mb} HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -L {params.regions} --native-pair-hmm-threads {threads} {params.extra}
		"""

rule GenomicsDB:
	"""
	Here GenomicsDB replace CombineGVCFs because it was taking too much time to run. GenomicsDB create a database and  store all variations present in gvcf file. Later sample will be genotyped
	"""
	input:
		expand("results/03_snv_indels_calling/{s.sample}.{{regions}}.g.vcf",s=SAMPLE.itertuples()),
	output:
		dir=temp(directory("results/03_snv_indels_calling/genomicDB_{regions}")),
		json="results/03_snv_indels_calling/genomicDB_{regions}/callset.json",
		vcf="results/03_snv_indels_calling/genomicDB_{regions}/vcfheader.vcf",
		vidmap="results/03_snv_indels_calling/genomicDB_{regions}/vidmap.json"
	params:
		region=get_region,
		#region=get_region_genomicdb,
		extra=config["params_genomics_db"]
	conda:
		"../envs/gatk4.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	shell:
		"""
		if [ -d {output.dir} ]
		then
			rm -R {output.dir}
		fi
		variant=""
		for file in {input}
		do
			variant=$variant" -V "$file
		done
		gatk --java-options -Xmx{resources.mem_mb} GenomicsDBImport $variant --genomicsdb-workspace-path {output.dir} --intervals {params.region} {params.extra}
		"""

#def get_input_genotype_gvcf(wildcards):
#	if bool_region == True:
#		return expand("results/03_snv_indels_calling/genomicDB_regions/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"])
#	else:
#		return expand("results/03_snv_indels_calling/genomicDB_{regions}/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"],regions=LIST_REGION_HAPLOTYPECALLER),
#

rule Genotype_gvcf:
	"""
	Genotype samples using gatk GenotypeGVCFs
	"""
	input:
		ref="results/genome/"+REF_NAME+".fasta",
		dir="results/03_snv_indels_calling/genomicDB_{regions}/",
		#db=expand("results/03_snv_indels_calling/genomicDB/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"]),
		db=expand("results/03_snv_indels_calling/genomicDB_{regions}/{file}",file=["callset.json","vcfheader.vcf","vidmap.json"],regions=LIST_REGION_HAPLOTYPECALLER),
	output:
		vcf=temp("results/03_snv_indels_calling/variants.{regions}.vcf"),
		idx=temp("results/03_snv_indels_calling/variants.{regions}.vcf.idx"),
	conda:
		"../envs/gatk4.yaml"
	resources:
		mem_mb=config["ram"]
	threads: get_thread
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
#gatk GenotypeGVCFs --reference {input.ref} -O {output} -V $list_db
#gatk GenotypeGVCFs --reference {input.ref} -V gendb://{input.dir} -O {output}

rule gatherVCF:
	input:
		vcf=expand("results/03_snv_indels_calling/variants.{regions}.vcf",regions=LIST_REGION_HAPLOTYPECALLER),
		idx=expand("results/03_snv_indels_calling/variants.{regions}.vcf.idx",regions=LIST_REGION_HAPLOTYPECALLER)
	output:
		temp("results/03_snv_indels_calling/variants.merged.vcf")
	conda:
		"../envs/picard.yaml"
	resources:
		mem_mb=config["ram"]
	threads: get_thread
	params:
		extra=config["params_gathervcf"]
	shell:
		"""
		list=""
		for file in {input.vcf}
		do
			list=$list" I="$file
		done
		picard -Xmx{resources.mem_mb} GatherVcfs $list O={output} {params.extra}
		"""
