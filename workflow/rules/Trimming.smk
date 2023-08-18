rule TrimReads:
	input:
		read1=lambda wildcards: f"{config['samples'][wildcards.sample]}_1.fastq",
		read2=lambda wildcards: f"{config['samples'][wildcards.sample]}_2.fastq"
	output:
		read1="{sample}/data/{sample}_1.trimmed.fastq",
		read2="{sample}/data/{sample}_2.trimmed.fastq",
		read1_un="{sample}/data/{sample}_U1.trimmed.fastq",
		read2_un="{sample}/data/{sample}_U2.trimmed.fastq"
	params:
		adapters=config["adapter_trimming"]
	resources:
		mem_mb=5000,
		time="0-02:00:00"
	envmodules:
		"Trimmomatic"
	shell:
		"""
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 {input.read1} {input.read2} {output.read1} {output.read1_un} {output.read2} {output.read2_un} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36
		"""

rule TrimmingQC_fastqc:
	input:
		read="{sample}/data/{sample}_{read}.trimmed.fastq"
	output:
		read="{sample}/qc/trimming/{sample}_{read}.trimmed.fastqc.html" 
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
		mkdir -p {wildcards.sample}/qc/trimming
		fastqc -o {wildcards.sample}/qc/trimming {input.read} 
		mv {params.read} {output.read}
		"""

rule TrimmingQC_multiqc:
	input:
		expand("{{sample}}/qc/trimming/{{sample}}_{read}.trimmed.fastqc.html", read=["1","2"])
	output:
		"{sample}/qc/trimming/trimming_report_{sample}.html"
	resources:
		mem_mb=5000,
		time="0-00:10:00"
	envmodules:
		"MultiQC"
	shell:
		"""
		multiqc -n {output} --no-data-dir {wildcards.sample}/qc/trimming
		"""