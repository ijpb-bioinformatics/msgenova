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

def get_thread(wildcards):
	return int(config["cpu"])

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
