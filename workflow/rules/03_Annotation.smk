rule GeSeqAutomation_AnnotationPreStandardization:
	input:
		"{sample}/02_assembly/{sample}.assembled.fasta",
	output:
		"{sample}/03_annotation/{sample}.annotated.gb"
	log:
		"{sample}/03_annotation/logs/{sample}_03_1GeSeqAutomation_AnnotationPreStandardization.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardize"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardize"]["time"],
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

rule StandardizationAnnotation:
	input:
		"{sample}/03_annotation/{sample}.annotated.gb",
	output:
		"{sample}/03_annotation/{sample}.standardized.fasta"
	log:
		"{sample}/03_annotation/logs/{sample}_03_2StandardizedAnnotation.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardize"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardize"]["time"]
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/StandardizationAnnotation.py"

rule GeSeqAutomation_AnnotationPostStandardization:
	input:
		"{sample}/03_annotation/{sample}.standardized.fasta",
	output:
		"{sample}/03_annotation/{sample}.standardized.gb"
	log:
		"{sample}/03_annotation/logs/{sample}_03_3GeSeqAutomation_AnnotationPostStandardization.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardize"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardize"]["time"],
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

def AnnotationQualityControlInputFunc(wildcards):
	if (not config["Standardization"]):
		return ["{sample}/02_assembly/{sample}.assembled.fasta","{sample}/03_annotation/{sample}.annotated.gb"]
	else:
		return ["{sample}/03_annotation/{sample}.standardized.fasta","{sample}/03_annotation/{sample}.standardized.gb"]


rule AnnotationQualityControl:
	input:
		AnnotationQualityControlInputFunc
	output:
		"{sample}/03_annotation/{sample}.standardized.filt.gb" if config["Standardization"] else "{sample}/03_annotation/{sample}.annotated.filt.gb"
	params:
		config["plant_genes_qc_file"] if os.path.exists(config["Standardize"]["plant_genes_qc_file"]) else os.path.join(workflow.basedir, config["Standardize"]["plant_genes_qc_file"]),
		alt_start=config["Standardize"]["alt_start"], 
		annotator=config["Standardize"]["annotator"]

	log:
		#"{sample}/03_annotation/logs/{sample}_03_4AnnotationQualityControl.log",
		"{sample}/03_annotation/{sample}.AnnotationCheckfile.log" if config["Standardization"] else "{sample}/03_annotation/logs/{sample}_03_4AnnotationQualityControl.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardize"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardize"]["time"]
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/AnnotationQualityControl.py"
