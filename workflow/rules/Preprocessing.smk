rule StartQC_fastqc:
	input:
		#read ="{sample}/data/{sample}_{read}.fastq"
		read = lambda wildcards: f"{config['samples'][wildcards.sample]}_{wildcards.read}.fastq"
	output:
		read="{sample}/qc/start/{sample}_{read}.fastqc.html"
	params:
		read = lambda wildcards, output: output.read.replace('.fastqc', '_fastqc')
	resources:
		mem_mb=5000,
		time="0-00:30:00"
	envmodules:
		"FastQC/0.11.9-Java-11"
	shell:
		"""
		mkdir -p {wildcards.sample}/qc
		mkdir -p {wildcards.sample}/qc/start
		fastqc -o {wildcards.sample}/qc/start {input.read} 
		mv {params.read} {output.read}
		"""

rule StartQC_multiqc:
	input:
		expand("{{sample}}/qc/start/{{sample}}_{read}.fastqc.html", read=["1","2"]),
	output:
		"{sample}/qc/start/StartReportQC_{sample}.html" 
	resources:
		mem_mb=5000,
		time="0-00:10:00"
	envmodules:
		"MultiQC"
	shell:
		"""
		multiqc -n {output} --no-data-dir {wildcards.sample}/qc/start
		"""