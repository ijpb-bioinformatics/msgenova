def get_fq1(wildcards):
	return SAMPLE.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2(wildcards):
	return SAMPLE.loc[(wildcards.sample), ["fq2"]].dropna()

#create needed directories
def create_logsdir(path):
	try:
		os.makedirs(path,exist_ok=True)
	except OSError:
		logger.error("Creation of the directory %s failed" % path)
	else:
		logger.info("Successfully created the directory %s " % path)

#define function to get threads from config file or set defaults value
def get_thread(rule):
	List_of_rules=["sort_sam","mark_duplicate","build_snpeff_db","annotate_variant","deal_with_vector","bwa_index_reference","create_dict_reference","samtools_index_reference","create_input_pindel","merge_vcf_files","annotate_pindel_output","multiqc","trimmed_fqc2","trimmed_fqc1","concatenate_log_trimmomatic","flagstat","samtools_coverage_by_regions","samtools_coverage","samtools_coverage_final","samtools_depth","copy_config","copy_region","copy_sample_file","extract_flagstat","report","HaplotypeCaller","GenomicsDB","Genotype_gvcf","gatherVCF","cut_vector_file","copy_tdna","index_vector","prepare_reference_tdnascan","clean_and_delete_tdnascan"]
	if rule in List_of_rules:
		return int("1")
	else:
		try:
			value=config["cpu"]
			return int(config["cpu"])
		except:
			return int("16")	

#return int(config["cpu"])

def get_mem(rule):
	List_of_rules=["build_snpeff_db","deal_with_vector","create_input_pindel","merge_vcf_files","multiqc","trimming","trimmed_fqc2","trimmed_fqc1","concatenate_log_trimmomatic","flagstat","samtools_coverage_by_regions","samtools_coverage","samtools_coverage_final","samtools_depth","copy_config","copy_region","copy_sample_file","extract_flagstat","report","cut_vector_file","copy_tdna","index_vector","prepare_reference_tdnascan","clean_and_delete_tdnascan"]
	if rule in List_of_rules:
		return str("2G")
	else:
		if get_vector(config,"mem")[0] == "TRUE":
			return str(get_vector(config,"mem")[1])
		else:
			return "60G"

def get_size_insert(wildcards):
	try:
		value=config["insert_size"]
		return int(config["insert_size"])
	except:
		return int("500")

def transform_project_name(name):
	new_name=re.sub(r"[.,;\n\t\s]*","",name)
	return new_name

def get_DP_min(wildcards):
	if get_vector(config,"DP.min")[0] == "TRUE":
		return str(get_vector(config,"DP.min")[1])
	else:
		return "5"

def get_AR_min(wildcards):
	if get_vector(config,"AR.min")[0] == "TRUE":
		return str(get_vector(config,"AR.min")[1])
	else:
		return "0.2"

def get_AR(wildcards):
	if get_vector(config,"AR")[0] == "TRUE":
		return str("TRUE")
	else
		return str("FALSE")

def get_config_key(cfg,key,bool_var,dict):
	try:
		value=cfg[key]
		print("regions file exist")
		bool_var= True
		#Define Regions for analysis
		LIST_CHR=[]
		dict={}
		dict[str(0)]=config["regions"]
		LIST_CHR.append(0)
		return LIST_CHR,dict
	except:
		print("no regions file - Read reference")
		bool_var= False
		default_list=read_reference(config["reference"])
		return default_list


def define_reference(cfg,key):
	try:
		value=cfg[key]
		print("there os vector in the config file")
		VECTOR_NAME=Path(config["vector"]).stem
		REFERENCE=REF_NAME+"_"+VECTOR_NAME
		return REFERENCE
	except:
		print("No vector in config")
		REFERENCE=REF_NAME
		return REFERENCE

def get_vector(cfg,key):
	try:
		value=cfg[key]
		bool="TRUE"
		return (bool,value)
	except:
		bool="FALSE"
		return (bool,None)


def read_reference(reference):
	LIST_CHR_REFERENCE=[]
	dict={}
	i=0
	# try to get the chromosome from fasta
	with open(reference,"r") as file:
	#with open(config["reference"],"r") as file:
		for line in file:
			if re.search(">",line):
				#if ">" is present then cut by space, and get first sequence
				split_line=line.split(" ")
				result=re.sub(">","",split_line[0])
				LIST_CHR_REFERENCE.append(i)
				dict[str(i)]=re.sub(r"[\n\t\s]*", "", result)
				i=i+1
	return LIST_CHR_REFERENCE,dict


def get_regions(wildcards):
	return wilcards.chr+":"+wilcards.beg+"-"+wilcards.end
