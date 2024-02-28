def annotationSubmissionNCBIInphutFunc(wildcards):
	output = []
	#for i in config["samples"]:
	if wildcards.sample in config["samples_metadata"]:
		output += [config["samples_metadata"][wildcards.sample]]
	if (not config["Standardization"]):
		return output + ["{sample}/03_annotation/{sample}.annotated.gb","{sample}/03_annotation/{sample}.assembled.fasta"]
	else:
		return output + ["{sample}/03_annotation/{sample}.standardardized.filt.gb","{sample}/03_annotation/{sample}.standardardized.fasta"]


rule annotationSubmissionNCBI:
	input:
		gb=annotationSubmissionNCBIInphutFunc
	output:
		"{sample}/04_submission/{sample}.submission.gb"
	resources:
		mem_mb=5000,
		time="0-0:20:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/AnnotationSubmission.py"