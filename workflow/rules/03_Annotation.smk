rule AnnotationPreStandardization:
	input:
		"{sample}/02_assembly/{sample}.assembled.fasta"
	output:
		"{sample}/03_annotation/{sample}.annotated.gb"
	resources:
		mem_mb=2000,
		time="0-00:30:00",
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

rule StandardizationAnnotation:
	input:
		"{sample}/03_annotation/{sample}.annotated.gb"
	output:
		"{sample}/03_annotation/{sample}.standardardized.fasta"
	resources:
		mem_mb=1000,
		time="0-0:20:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/StandardizationAnnotation.py"

rule AnnotationPostStandardization:
	input:
		"{sample}/03_annotation/{sample}.standardardized.fasta"
	output:
		"{sample}/03_annotation/{sample}.standardardized.gb"
	resources:
		mem_mb=2000,
		time="0-00:30:00",
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

if (not config["Standardization"]):
	rule AnnotationQualityControl:
		input:
			"{sample}/02_assembly/{sample}.assembled.fasta",
			"{sample}/03_annotation/{sample}.annotated.gb",
		output:
			"{sample}/03_annotation/{sample}.annotated.filt.gb"
		params:
			config["plant_genes_file"]
		log:
			"{sample}/03_annotation/{sample}.annotated.warnings.log"
		resources:
			mem_mb=5000,
			time="0-0:20:00"
		conda:
			"../envs/Annotation.yaml"
		script:
			"../scripts/AnnotationQualityControl.py"

else:
	rule AnnotationQualityControl:
		input:
			"{sample}/03_annotation/{sample}.standardardized.fasta",
			"{sample}/03_annotation/{sample}.standardardized.gb"
		output:
			"{sample}/03_annotation/{sample}.standardardized.filt.gb"
		params:
			config["plant_genes_file"]
		log:
			"{sample}/03_annotation/{sample}.standardardized.warnings.log"
		resources:
			mem_mb=5000,
			time="0-0:30:00"
		conda:
			"../envs/Annotation.yaml"
		script:
			"../scripts/AnnotationQualityControl.py"
