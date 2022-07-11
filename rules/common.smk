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
