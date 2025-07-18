def buildReferenceAnnotationBackmappingInputFunc(wildcards):
	if (config["Standardization"]):
		return "{sample}/03_annotation/{sample}.standardized.fasta"
	else:
		return "{sample}/02_assembly/{sample}.assembled.fasta"

rule BuildReferenceAnnotationBackmapping:
	input:
		ref=buildReferenceAnnotationBackmappingInputFunc
	output:
		"{sample}/04_backmapping/mapping_ref/{sample}.backmapping.1.bt2"
	log:
		"{sample}/04_backmapping/qc/{sample}.bowtie2.reference.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*2/1000000),config["Backmapping"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Backmapping"]["time"]
	threads:8
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		bowtie2-build --threads {threads} {input.ref} {wildcards.sample}/04_backmapping/mapping_ref/{wildcards.sample}.backmapping 2> {log}
		"""

rule Backmapping:
	input:
		read1="{sample}/01_data/{sample}_1.trimmed.fastq.gz",
		read2="{sample}/01_data/{sample}_2.trimmed.fastq.gz",
		ref="{sample}/04_backmapping/mapping_ref/{sample}.backmapping.1.bt2"
	output:
		"{sample}/04_backmapping/mapping/{sample}.backmapping.sam"
	log:
		"{sample}/04_backmapping/qc/{sample}.bowtie2.backmapping.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*3/1000000),config["Backmapping"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Backmapping"]["time"],
	threads:8
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		bowtie2 --threads {threads} -x {wildcards.sample}/04_backmapping/mapping_ref/{wildcards.sample}.backmapping -1 {input.read1} -2 {input.read2} -S {output} 2> {log}
		"""

rule BackmappingQC:
	input:
		"{sample}/04_backmapping/mapping/{sample}.backmapping.sam"
	output:
		sbam="{sample}/04_backmapping/{sample}.backmapping.sorted.bam",
		qc="{sample}/04_backmapping/{sample}.backmapping.readcoverage.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*1/10000000),config["Backmapping"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Backmapping"]["time"]
	conda:
		"../envs/Samtools.yaml"
	shell:
		"""
		samtools view -bS {input} | samtools sort > {output.sbam}
		samtools coverage {output.sbam} > {output.qc}
		samtools coverage -A -w 64 {output.sbam} >> {output.qc}
		"""

rule BackmappingGenomeCoverage:
	input:
		sbam="{sample}/04_backmapping/{sample}.backmapping.sorted.bam"
	output:
		cov="{sample}/04_backmapping/qc/{sample}.backmapping.coverage.genomecov",
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*2/1000000),config["Backmapping"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Backmapping"]["time"]
	conda:
		"../envs/Bedtools.yaml"
	shell:
		"""
		bedtools genomecov -d -ibam {input.sbam} > {output.cov}
		"""


rule BackmappingVisualization:
	input:
		"{sample}/04_backmapping/qc/{sample}.backmapping.coverage.genomecov",
	output:
		vis="{sample}/04_backmapping/{sample}.backmapping.visualization.html"
	log:
		"{sample}/04_backmapping/qc/{sample}.backmapping.low_coverage.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[0]).st_size*2/1000000),config["Backmapping"]["mem_mb"]),
		time=lambda wildcards, attempt: config["Backmapping"]["time"]
	conda:
		"../envs/Annotation.yaml"
	script:
		"../scripts/BackmappingVisualization.py"
