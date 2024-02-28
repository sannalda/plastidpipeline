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

	rule AnnotationQualityCheck:
		input:
			"{sample}/03_annotation/{sample}.annotated.filt.gb"
		output:
			"{sample}/03_annotation/{sample}.unpause.txt"
		resources:
			chrome=1
		localrule: True
		shell:
			"""
			echo "Hello, please introduce yourself."

			echo -n "Your name: "
			read -r name
			"""

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

	rule AnnotationQualityCheck:
		input:
			sfilt="{sample}/03_annotation/{sample}.standardardized.filt.gb",
			log="{sample}/03_annotation/{sample}.standardardized.warnings.log"
		output:
			"{sample}/03_annotation/{sample}.unpause.txt"
		resources:
			chrome=1
		localrule: True
		shell:
			"""
			echo "Annoation finished for {wildcards.sample}. Please check {input.log} and make any necessary changes to {input.sfilt} before proceeding with the pipeline. 
			To proceed, type 'y'. Else, type anything (i.e. 'n') and enter to exit the pipeline. This diaglogue will repear again when you rerun the pipeline if you exit."

			echo -n "Response (y or n): "
			read -r response

			if [ "$response" = "y" ]; then 
				touch {output}
			else 
				echo "Exiting the pipeline in 60 seconds..." 
			fi
			"""
