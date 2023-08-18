def buildReferenceAnnotationBackmappingInputFunc(wildcards):
	if (config["Standardization"]):
		return rules.StandardizationAnnotation.output
	else:
		return rules.Assembly.output.assembly

rule BuildReferenceAnnotationBackmapping:
	input:
		buildReferenceAnnotationBackmappingInputFunc
	output:
		"{sample}/backmapping/{sample}.backmapping.1.bt2"
	log:
		"{sample}/qc/backmapping/{sample}.bowtie2.reference.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	threads:4
	envmodules:
		"Bowtie2"
	shell:
		"""
		bowtie2-build --threads {threads} {input} {wildcards.sample}/backmapping/{wildcards.sample}.backmapping 2> {log}
		"""

rule Backmapping:
	input:
		ref="{sample}/backmapping/{sample}.backmapping.1.bt2",
		read1="{sample}/data/{sample}_1.filt.fastq",
		read2="{sample}/data/{sample}_2.filt.fastq"
	output:
		"{sample}/backmapping/{sample}.backmapping.sam"
	log:
		"{sample}/qc/backmapping/{sample}.bowtie2.backmapping.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	threads:8
	envmodules:
		"Bowtie2"
	shell:
		"""
		bowtie2 --threads {threads} -x {wildcards.sample}/backmapping/{wildcards.sample}.backmapping -1 {input.read1} -2 {input.read2} -S {output} 2> {log}
		"""

rule BackmappingQC:
	input:
		rules.Backmapping.output
	output:
		bam="{sample}/backmapping/{sample}.backmapping.bam",
		sbam="{sample}/backmapping/{sample}.backmapping.sorted.bam",
		qc="{sample}/qc/backmapping/{sample}.backmapping.readcoverage.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	envmodules:
		"SAMtools"
	shell:
		"""
		samtools view -bS {input} > {output.bam}
		samtools sort {output.bam} > {output.sbam}
		samtools coverage {output.sbam} > {output.qc}
		samtools coverage -A -w 64 {output.sbam} >> {output.qc}
		"""