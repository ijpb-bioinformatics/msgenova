
rule copy_config:
	"""
	Copy config advanced in 00_logs
	"""
	output:
		"results/00_logs/config_advanced"
	params:
		config=config["repo_script"]
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"""
		cp {params.config}/config_advanced {output}
		"""

rule report:
	"""
	Create final report comporting all analysis results
	"""
	input:
		multiqc="results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html",
		file_depth=expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),
		#flag="results/03_mapping/flagstat/Extract_data.flagstat",
		flag=expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),
		sample_file=config["sample"],
		file_coverage=expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		#expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		reference=config["reference"],
		qc_trimmo="results/01_sequence_qc/log/trimmomatic.log",
		haplotype_caller_html="results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html",
		haplotype_caller_vcf="results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf",
		config="results/00_logs/config_advanced",
		insertion=aggregate_vector,
		#pindel=get_input_gather_pindel
	output:
		"results/06_report/msgenova_report.html"
	params:
		workdir=config["repo_script"],
		DP=config["DP.min"],
		AR=config["AR.min"],
		vector=get_vector(config,"vector")[1]
	conda:
		"../envs/R.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"Rscript -e \"rmarkdown::render('{params.workdir}/script/ms_report.Rmd', output_file = '{params.workdir}/{output}', params = list(result_dir='{params.workdir}/results/', sample='{input.sample_file}', DP.min='{params.DP}', AR.min='{params.AR}', multiqc='{params.workdir}/{input.multiqc}', vector='{params.vector}', trimmo='{params.workdir}/{input.qc_trimmo}',flagstat='{params.workdir}/{input.flag}',haplotype_caller_vcf='{params.workdir}/{input.haplotype_caller_vcf}'))\""

