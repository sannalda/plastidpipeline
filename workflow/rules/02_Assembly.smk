rule AssemblyConfig:
	output:
		expand("%s/SeedDatabase/{organelle_dbs}.fasta" %config["organelle_database_folder"], organelle_dbs=config["OrganelleDatabasesDownload"].split(","))
	resources:
		mem_mb=20000,
		time="0-02:00:00"
	envmodules:
		"Bowtie2",
		"BLAST+"
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		mkdir -p {config[organelle_database_folder]}
		get_organelle_config.py -a {config[OrganelleDatabasesDownload]} --config-dir {config[organelle_database_folder]}
		"""

rule Assembly:
	input:
		expand("%s/SeedDatabase/{organelle_dbs}.fasta" %config["organelle_database_folder"], organelle_dbs=config["OrganelleDatabasesAssembly"].split(",")),
		read1="{sample}/01_data/{sample}_1.trimmed.sorted.fastq.gz",
		read2="{sample}/01_data/{sample}_2.trimmed.sorted.fastq.gz"
	output:
		assembly="{sample}/02_assembly/{sample}.assembled.fasta"
	params:
		assembly_pattern=".complete.graph1.1.path_sequence.fasta"
	resources:
		mem_mb=20000,
		time="0-6:00:00"
	threads: 8
	envmodules:
		"Bowtie2",
		"BLAST+",
		"SPAdes"
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		get_organelle_from_reads.py -1 {input.read1} -2 {input.read2} \\
			-o {wildcards.sample}/02_assembly -R {config[MaxRounds]} -k {config[KmerSpades]} -P {config[PreGrouped]} \\
			-F {config[OrganelleDatabasesAssembly]} --config-dir {config[organelle_database_folder]}\\
			-t {threads} --prefix {wildcards.sample} --overwrite
		scp {wildcards.sample}/02_assembly/*{params.assembly_pattern} {output.assembly}
		sed -i '1s/.*/>>{wildcards.sample}/' {output.assembly}
		"""
