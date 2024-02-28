def buildReferenceAnnotationBackmappingInputFunc(wildcards):
	if (config["Standardization"]):
		return "{sample}/03_annotation/{sample}.standardardized.fasta"
	else:
		return "{sample}/02_assembly/{sample}.assembled.fasta"

rule BuildReferenceAnnotationBackmapping:
	input:
		buildReferenceAnnotationBackmappingInputFunc
	output:
		"{sample}/04_backmapping/mapping/{sample}.backmapping.1.bt2"
	log:
		"{sample}/04_backmapping/qc/{sample}.bowtie2.reference.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	threads:4
	envmodules:
		"Bowtie2"
	shell:
		"""
		bowtie2-build --threads {threads} {input} {wildcards.sample}/04_backmapping/mapping/{wildcards.sample}.backmapping 2> {log}
		"""

rule Backmapping:
	input:
		ref="{sample}/04_backmapping/mapping/{sample}.backmapping.1.bt2",
		read1="{sample}/01_data/{sample}_1.filt.fastq",
		read2="{sample}/01_data/{sample}_2.filt.fastq"
	output:
		"{sample}/04_backmapping/mapping{sample}.backmapping.sam"
	log:
		"{sample}/04_backmapping/qc/{sample}.bowtie2.backmapping.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	threads:8
	envmodules:
		"Bowtie2"
	shell:
		"""
		bowtie2 --threads {threads} -x {wildcards.sample}/04_backmapping/mapping/{wildcards.sample}.backmapping -1 {input.read1} -2 {input.read2} -S {output} 2> {log}
		"""

rule BackmappingQC:
	input:
		"{sample}/04_backmapping/mapping/{sample}.backmapping.sam"
	output:
		bam="{sample}/04_backmapping/mapping/{sample}.backmapping.bam",
		sbam="{sample}/04_backmapping/{sample}.backmapping.sorted.bam",
		qc="{sample}/04_backmapping/{sample}.backmapping.readcoverage.log"
	resources:
		mem_mb=10000,
		time="0-01:00:00"
	conda:
		"../envs/Samtools.yaml"
	shell:
		"""
		samtools view -bS {input} > {output.bam}
		samtools sort {output.bam} > {output.sbam}
		samtools coverage {output.sbam} > {output.qc}
		samtools coverage -A -w 64 {output.sbam} >> {output.qc}
		"""

rule BackmappingGenomeCov:
	input:
		sbam="{sample}/04_backmapping/{sample}.backmapping.sorted.bam"
	output:
		cov="{sample}/04_backmapping/qc/{sample}.backmapping.coverage.genomecov",
	resources:
		mem_mb=10000,
		time="0-00:30:00"
	conda:
		"../envs/Bedtools.yaml"
	shell:
		"""
		bedtools genomecov -d -ibam {input.sbam} > {output.cov}
		"""


rule BackmappingVisualization:
	input:
		"{sample}/04_backmapping/qc/{sample}.backmapping.coverage.genomecov"
	output:
		vis="{sample}/04_backmapping/{sample}.backmapping.visualization.html"
	log:
		"{sample}/04_backmapping/qc/{sample}.backmapping.low_coverage.log"
	resources:
		mem_mb=1000,
		time="0-00:15:00"
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/BackmappingVisualization.py"
