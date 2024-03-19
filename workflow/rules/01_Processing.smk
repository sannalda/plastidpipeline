rule TrimReads:
	input:
		read1=lambda wildcards: f"{config['samples'][wildcards.sample]}_1"+end,
		read2=lambda wildcards: f"{config['samples'][wildcards.sample]}_2"+end
	output:
		read1="{sample}/01_data/{sample}_1.trimmed"+end,
		read2="{sample}/01_data/{sample}_2.trimmed"+end,
		qc="{sample}/01_data/{sample}.trimmed.qc.html",
		json="{sample}/01_data/{sample}.trimmed.qc.json"
	params:
		adapters=config["adapter_trimming_file"]
	resources:
		mem_mb=20000,
		time="0-02:00:00"
	threads: 8
	conda:
		"../envs/Fastp.yaml"
	shell:
		"""
		fastp -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} -h {output.qc} -j {output.json} --detect_adapter_for_pe --adapter_fasta {params.adapters} -w {threads}
		"""


def sortReadsInputFunc(wildcards):
	if (config["Trimming"]):
		return "{sample}/01_data/{sample}_{read}.trimmed"+end
	else:
		return f"{config['samples'][wildcards.sample]}_{wildcards.read}{end}"

rule SortReads:
	input:
		read=sortReadsInputFunc			
	output:
		read="{sample}/01_data/{sample}_{read}.trimmed.sorted.fastq.gz"
	resources:
		mem_mb=20000,
		time="0-03:00:00"
	envmodules:
		"bioawk"
	shell:
		"""
		bioawk -c fastx '{{print}}' {input} | sort | awk -F'\\t' '{{print "@"$1;print $2;print "+"$1;print $3}}' | gzip > {output.read}
		"""
