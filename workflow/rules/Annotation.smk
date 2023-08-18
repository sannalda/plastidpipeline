rule AnnotationPreStandardization:
	input:
		rules.Assembly.output.assembly
	output:
		"{sample}/annotation/{sample}.original.gb"
	conda:
		"../envs/test.yaml"
	script:
		"../scripts/GeSeqAutomation.py"

rule StandardizationAnnotation:
	input:
		rules.AnnotationPreStandardization.output
	output:
		"{sample}/annotation/{sample}.standardardized.fasta"
	resources:
		mem_mb=1000,
		time="0-0:30:00"
	conda:
		"../envs/test.yaml"
	script:
		"../scripts/StandardizationAnnotation.py"

rule AnnotationPostStandardization:
	input:
		rules.StandardizationAnnotation.output
	output:
		"{sample}/annotation/{sample}.standardardized.gb"
	conda:
		"../envs/test.yaml"
	script:
		"../scripts/GeSeqAutomation.py"