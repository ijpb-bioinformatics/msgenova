rule report:
	"""
	Create final report comporting all analysis results
	"""
	input:
		multiqc="results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html",
		#file_depth=expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),
		#flag="results/03_mapping/flagstat/Extract_data.flagstat",
		flag=expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),
		sample_file=config["sample"],
		file_coverage=expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		reference=config["reference"],
		haplotype_caller_html="results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html",
		haplotype_caller_vcf="results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf",
	output:
		"results/06_report/msgenova_report.html"
	params:
		workdir=config["repo_script"],
		DP=4,
		AR=0.2,
		region=dict_chr.values()
	conda:
		"../envs/R.yaml"
	threads: get_thread
	resources:
		mem_mb=config["ram"]
	shell:
		"Rscript -e \"rmarkdown::render('{params.workdir}/script/ms_report.Rmd', output_file = '{params.workdir}/{output}', params = list(result_dir='{params.workdir}/results/', sample='{input.sample_file}',region='{params.region}', DP.min='{params.DP}', AR.min='{params.AR}', multiqc='{params.workdir}/{input.multiqc}'))\""

