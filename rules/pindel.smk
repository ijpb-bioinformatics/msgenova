rule create_input_pindel:
	input:
		expand("results/02_mapping/bam/{s.sample}.bam", s=SAMPLE.itertuples()),
	output:
		"results/04_pindel/input_pindel.txt"
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
		for file in {input}
		do
			nom=`basename $file`
			echo -e "{params.wd}/$file\t{params.ins_size}\t$nom" >> {output}
		done
		"""

def get_regions_pindel(wildcards):
	return wildcards.region_pindel
#	#if there is a region file, construct region, else get chromosome
#	if (get_vector(config,"regions")[0] == "TRUE") :
#		#then construct regions
#		return widcards.chr+":"+wildcards.beg+"-"+wildcards.end
#	else:
#		return wildcards.region_pindel
#		#else there are no regions file, so get chromosome list
#		#x=dict_chr.get(str(wildcards.region_pindel))
#		#return x


rule run_pindel:
	input:
		config_pindel="results/04_pindel/input_pindel.txt",
		files=expand("results/02_mapping/bam/{s.sample}.bam", s=SAMPLE.itertuples()),
		reference="results/genome/"+REFERENCE+".fasta",
		index_bam=expand("results/02_mapping/bam/{s.sample}.bam.bai", s=SAMPLE.itertuples()),
		fai_ref="results/genome/"+REFERENCE+".fasta.fai"
	output:
		expand("results/04_pindel/{{region_pindel}}_{ext_pindel}",ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),
	threads: get_thread
	resources: mem_mb=get_mem
	params:
		wd=WORKDIR,
		reg=get_regions_pindel
	conda:
		"../envs/pindel.yaml"
	shell:
		"""
		pindel -T {threads} -f {input.reference} -i {input.config_pindel} -c {params.reg} -o {params.wd}/results/04_pindel/{params.reg}
		"""


def get_input_gather_pindel(wilkdcards):
	list=[]
	if (get_vector(config,"regions")[0] == "TRUE") :
		region_p=[]
		for elm in REGIONS.itertuples():
			reg=str(elm.chr+":"+elm.beg+"-"+elm.end)
			region_p.append(reg)
	else:
		region_p=list(dict_chr.values())
	list.extend(expand("results/04_pindel/{region_pindel}_{ext_pindel}",region_pindel=region_p,ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),)
	print(list)
	return list

rule gather_pindel:
	input:
		get_input_gather_pindel
		#expand("results/04_pindel/{region_pindel}_{ext_pindel}",region_pindel=list(dict_chr.values()),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),
	output:
		"aggregate_pindel.txt"
	threads: get_thread
	resources: mem_mb=get_mem
	shell:
		"""
		for file in {input}
		do
			echo $file >> {output}
		done
		"""
