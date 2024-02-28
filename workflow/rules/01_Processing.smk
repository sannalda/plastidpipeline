def sortReadsInputFunc(wildcards):
	if (config["Trimming"]):
		return "{sample}/01_data/{sample}_{read}.trimmed.fastq"
	else:
		return f"{config['samples'][wildcards.sample]}_{wildcards.read}.fastq"

rule SortReads:
	input:
		read=sortReadsInputFunc			
	output:
		read="{sample}/01_data/{sample}_{read}.sorted.fastq"
	resources:
		mem_mb=5000,
		time="0-01:00:00"
	envmodules:
		"bioawk"
	shell:
		"""
		bioawk -c fastx '{{print}}' {input} | sort | awk -F'\\t' '{{print "@"$1;print $2;print "+"$1;print $3}}' > {output.read}
		"""

rule FilteringReads:
	input:
		read=rules.SortReads.output.read
	output:
		read="{sample}/01_data/{sample}_{read}.filt.fastq"
	params:
		q=config["MinQualityScore"]
	resources:
		mem_mb=10000,
		time="0-03:00:00"
	envmodules:
		"FASTX-Toolkit"
	shell:
		"""
		fastq_quality_filter -q {config[MinQualityScore]} -i {input.read} -o {output.read}
		"""

rule FilteringQC_fastqc:
	input:
		read = rules.FilteringReads.output.read
	output:
		read = "{sample}/01_data/qc/filtering/{sample}_{read}.filt.fastqc.html"
	params:
		read = lambda wildcards, output: output.read.replace('.fastqc', '_fastqc'),
	resources:
		mem_mb=5000,
		time="0-00:30:00"
	envmodules:
		"FastQC/0.11.9-Java-11"
	shell:
		"""
		mkdir -p {wildcards.sample}/01_data/qc
		mkdir -p {wildcards.sample}/01_data/qc/filtering
		fastqc -o {wildcards.sample}/01_data/qc/filtering {input.read}
		mv {params.read} {output.read}
		"""

rule FilteringQC_multiqc:
	input:
		expand("{{sample}}/01_data/qc/filtering/{{sample}}_{read}.filt.fastqc.html", read=["1","2"]),
	output:
		"{sample}/01_data/qc/filtering_report_{sample}.html" 
	resources:
		mem_mb=5000,
		time="0-00:10:00"
	envmodules:
		"MultiQC"
	shell:
		"""
		multiqc -n {output} --no-data-dir {wildcards.sample}/01_data/qc/filtering
		"""
