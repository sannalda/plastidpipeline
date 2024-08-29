rule GeSeqAutomation_AnnotationPreStandardization:
	input:
		"{sample}/02_assembly/{sample}.assembled.fasta",
	output:
		"{sample}/03_annotation/{sample}.annotated.gb"
	log:
		"{sample}/logs/{sample}_03_1GeSeqAutomation_AnnotationPreStandardization.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardization"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardization"]["time"],
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
		"{sample}/03_annotation/{sample}.standardardized.fasta"
	log:
		"{sample}/logs/{sample}_03_2StandardizedAnnotation.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardization"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardization"]["time"]
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/StandardizationAnnotation.py"

rule GeSeqAutomation_AnnotationPostStandardization:
	input:
		"{sample}/03_annotation/{sample}.standardardized.fasta",
	output:
		"{sample}/03_annotation/{sample}.standardardized.gb"
	log:
		"{sample}/logs/{sample}_03_3GeSeqAutomation_AnnotationPostStandardization.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardization"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardization"]["time"],
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
		return ["{sample}/03_annotation/{sample}.standardardized.fasta","{sample}/03_annotation/{sample}.standardardized.gb"]


rule AnnotationQualityControl:
	input:
		AnnotationQualityControlInputFunc
	output:
		"{sample}/03_annotation/{sample}.standardardized.filt.gb" if config["Standardization"] else "{sample}/03_annotation/{sample}.annotated.filt.gb"
	params:
		config["plant_genes_file"] if os.path.exists(config["Standardization"]["plant_genes_qc_file"]) else os.path.join(workflow.basedir, config["Standardization"]["plant_genes_qc_file"])
	log:
		"{sample}/logs/{sample}_03_4AnnotationQualityControl.log",
		#"{sample}/03_annotation/{sample}.standardardized.warnings.log" if config["Standardization"] else "{sample}/03_annotation/{sample}.annotated.warnings.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Standardization"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Standardization"]["time"]
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/AnnotationQualityControl.py"
