rule run_delly:
	"""
	Run Delly sample by sample
	"""
	input:
		bam="results/02_mapping/bam/{sample}.bam",
		index_bam="results/02_mapping/bam/{sample}.bam.bai",
		reference="results/genome/"+REFERENCE+".fasta",
		fai_ref="results/genome/"+REFERENCE+".fasta.fai"
	output:
		bcf=temp("results/04_delly/{sample}.bcf"),
		csi=temp("results/04_delly/{sample}.bcf.csi")
	conda:
		"../envs/env_delly.yaml"
	shell:
		"""
		delly call -g {input.reference} -o {output.bcf}  {input.bam}
		"""

rule delly_bcf_2_vcf:
	"""
	convert bcf to vcf + index
	"""
	input:
		bcf="results/04_delly/{sample}.bcf"
	output:
		vcf= temp("results/04_delly/{sample}.vcf")
	conda:
		"../envs/env_alignment.yaml"
	shell:
		"""
		bcftools view {input.bcf} > {output.vcf}
		"""


rule vcf_2_vcfgz:
	"""
	compress vcf
	"""
	input:
		vcf="results/04_delly/{sample}.vcf"
	output:
		vcfgz="results/04_delly/{sample}.vcf.gz"
	conda:
		"../envs/env_alignment.yaml"
	shell:
		"""
		bgzip {input.vcf}
		"""
