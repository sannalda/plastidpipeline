def annotationSubmissionNCBIInphutFunc(wildcards):
	output = [config["metadata_file"]]
	if (not config["Standardization"]):
		return output + ["{sample}/annotation/{sample}.original.cleaned.gb","{sample}/annotation/{sample}.original.fasta"]
	else:
		return output + ["{sample}/annotation/{sample}.standardardized.cleaned.gb","{sample}/annotation/{sample}.standardardized.fasta"]


rule annotationSubmissionNCBI:
	input:
		gb=annotationSubmissionNCBIInphutFunc
	output:
		"{sample}/annotation/{sample}.submission.gb"
	resources:
		mem_mb=1000,
		time="0-0:20:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/AnnotationSubmission.py"