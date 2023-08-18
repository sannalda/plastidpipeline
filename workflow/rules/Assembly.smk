rule AssemblyConfig:
	output:
		expand("%s/SeedDatabase/{organelle_dbs}.fasta" %config["organelle_database_folder"], organelle_dbs=config["OrganelleDatabasesDownload"].split(","))
	resources:
		mem_mb=10000,
		time="0-01:00:00"
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
		read1="{sample}/data/{sample}_1.filt.fastq",
		read2="{sample}/data/{sample}_2.filt.fastq"
	output:
		assembly="{sample}/assembly/{sample}.original.fasta"
	params:
		assembly_pattern=".complete.graph1.1.path_sequence.fasta" # This is a janky way of doing it, but it works haha
	resources:
		mem_mb=20000,
		time="0-5:00:00"
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
			-o {wildcards.sample}/assembly -R {config[MaxRounds]} -k {config[KmerSpades]} -P {config[PreGrouped]} \\
			-F {config[OrganelleDatabasesAssembly]} --config-dir {config[organelle_database_folder]}\\
			-t {threads} --prefix {wildcards.sample} --overwrite
		scp {wildcards.sample}/assembly/*{params.assembly_pattern} {output.assembly}
		"""
