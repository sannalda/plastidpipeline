def AnnotationSubmissionNCBIInputFunc(wildcards):
	output = []
	#for i in config["samples"]:
	if wildcards.sample in config["samples_metadata"]:
		output += [config["samples_metadata"][wildcards.sample]]
	if (not config["Standardization"]):
		return output + ["{sample}/03_annotation/{sample}.annotated.gb","{sample}/03_annotation/{sample}.assembled.fasta"]
	else:
		return output + ["{sample}/03_annotation/{sample}.standardized.filt.gb","{sample}/03_annotation/{sample}.standardized.fasta"]


rule AnnotationMetadata:
	input:
		gb=AnnotationSubmissionNCBIInputFunc
	output:
		"{sample}/05_metadata/{sample}.final.gb"
	log:
		"{sample}/05_metadata/{sample}_05_1AnnotationMetadata.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*2/1000000),1000),
		time="0-0:30:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/AnnotationMetadata.py"
