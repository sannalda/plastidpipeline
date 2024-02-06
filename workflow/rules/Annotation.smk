rule AnnotationPreStandardization:
	input:
		"{sample}/assembly/{sample}.original.fasta"
	output:
		"{sample}/annotation/{sample}.original.gb"
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
		mem_mb=2000,
		time="0-00:30:00",
		chrome=1
	conda:
		"../envs/Annotation.yaml"
	localrule: True
	script:
		"../scripts/GeSeqAutomation.py"

if (not config["Standardization"]):
	rule annotationQualityCheckInputFunc:
		input:
			"{sample}/assembly/{sample}.original.fasta",
			"{sample}/annotation/{sample}.original.gb"
		output:
			"{sample}/annotation/{sample}.original.cleaned.gb",
			"{sample}/annotation/{sample}.original.incorrect.gb"
		log:
			"{sample}/annotation/{sample}.original.incorrect.log"
		resources:
			mem_mb=1000,
			time="0-0:20:00"
		conda:
			"../envs/Annotation.yaml"
		script:
			"../scripts/AnnotationQualityControl.py"
else:
	rule annotationQualityCheckInputFunc:
		input:
			"{sample}/annotation/{sample}.standardardized.fasta",
			"{sample}/annotation/{sample}.standardardized.gb"
		output:
			"{sample}/annotation/{sample}.standardardized.cleaned.gb",
			"{sample}/annotation/{sample}.standardardized.incorrect.gb"
		log:
			"{sample}/annotation/{sample}.standardardized.incorrect.log"
		resources:
			mem_mb=1000,
			time="0-0:20:00"
		conda:
			"../envs/Annotation.yaml"
		script:
			"../scripts/AnnotationQualityControl.py"