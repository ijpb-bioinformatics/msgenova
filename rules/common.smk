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
def get_thread(wildcards):
	try:
		value=config["cpu"]
		return int(config["cpu"])
	except:
		return int("16")	
#return int(config["cpu"])

def get_mem(wildcards):
	if get_vector(config,"ram")[0] == "TRUE":
		return str(get_vector(config,"ram")[1])
	else:
		return "30G"
#	try:
#		value=config["ram"]
#		return config["ram"]
#	except:
#		return "30G"
#

def transform_project_name(name):
	new_name=re.sub(r"[.,;\n\t\s]*","",name)
	return new_name

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
