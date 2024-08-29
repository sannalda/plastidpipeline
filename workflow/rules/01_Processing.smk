class InputDataFormat(Exception):
    pass

def DetermineFileExt(sample_name):
	if  os.path.exists("%s_1.fq.gz" %config["samples"][sample_name]):
		end=".fq.gz"
	elif os.path.exists("%s_1.fastq.gz" %config["samples"][sample_name]):
		end=".fastq.gz"
	elif  os.path.exists("%s_1.fastq" %config["samples"][sample_name]):
		end=".fastq"
	elif os.path.exists("%s_1.fq" %config["samples"][sample_name]):
		end=".fq"
	else:
		raise InputDataFormat
	return end

#def AdjustResources_TrimReads(wildcards, attempt):
    # Adjust resources based on the number of attempts
#    if attempt == 1:
#    	print("Trying at starting time limit 1 (00:30:00)...")
#        return {"time": "00:30:00"}
#    elif attempt == 2:
#    	print("Trying at increased time limit 2 (02:00:00)...")
#        return {"time": "02:00:00"}
#    elif attempt >= 3:
#    	print("Trying at increased time limit 3 (%s). If this doesn't work, try changing 'time' parameter in config file..." %config["time"])
#        return {"time": config["time"]}

rule TrimReads:
	input:
		read1=lambda wildcards: f"{config['samples'][wildcards.sample]}_1"+DetermineFileExt(wildcards.sample),
		read2=lambda wildcards: f"{config['samples'][wildcards.sample]}_2"+DetermineFileExt(wildcards.sample)
	output:
		read1="{sample}/01_data/{sample}_1.trimmed.fastq.gz",
		read2="{sample}/01_data/{sample}_2.trimmed.fastq.gz",
		qc="{sample}/01_data/{sample}.trimmed.qc.html",
		json="{sample}/01_data/{sample}.trimmed.qc.json"
	log:
		"{sample}/01_data/{sample}_01_TrimReads.log"
	params:
		adapters = config["Trimming"]["adapter_trimming_file"] if os.path.exists(config["Trimming"]["adapter_trimming_file"]) else os.path.join(workflow.basedir, config["Trimming"]["adapter_trimming_file"]),
		additional_params = config["Trimming"]["Fastp_AdditionalParameters"]
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Trimming"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Trimming"]["time"]
	threads: 16
	conda:
		"../envs/Fastp.yaml"
	shell:
		"""
		fastp -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} -h {output.qc} -j {output.json} --adapter_fasta {params.adapters} -w {threads} -V {params.additional_params} 2> {log}
		"""