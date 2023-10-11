rule AnnotationPreStandardization:
	input:
		"{sample}/assembly/{sample}.original.fasta"
	output:
		"{sample}/annotation/{sample}.original.gb"
	resources:
		mem_mb=1000,
		time="0-00:20:00",
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

rule StandardizationAnnotation:
	input:
		rules.AnnotationPreStandardization.output
	output:
		"{sample}/annotation/{sample}.standardardized.fasta"
	resources:
		mem_mb=1000,
		time="0-0:20:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/StandardizationAnnotation.py"

rule AnnotationPostStandardization:
	input:
		"{sample}/annotation/{sample}.standardardized.fasta"
	output:
		"{sample}/annotation/{sample}.standardardized.gb"
	resources:
		mem_mb=1000,
		time="0-00:20:00",
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"