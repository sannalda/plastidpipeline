rule AssemblyConfig:
	output:
		expand("%s/SeedDatabase/{organelle_dbs}.fasta" %config["Assembly"]["organelle_database_folder"], organelle_dbs=config["Assembly"]["OrganelleDatabasesDownload"].split(","))
	resources:
		mem_mb=5000,
		time="0-00:30:00"
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		mkdir -p {config[organelle_database_folder]}
		get_organelle_config.py -a {config[OrganelleDatabasesDownload]} --config-dir {config[organelle_database_folder]}
		"""

rule Assembly:
	input:
		expand("%s/SeedDatabase/{organelle_dbs}.fasta" %config["Assembly"]["organelle_database_folder"], organelle_dbs=config["Assembly"]["OrganelleDatabasesAssembly"].split(",")),
		read1="{sample}/01_data/{sample}_1.trimmed.fastq.gz",
		read2="{sample}/01_data/{sample}_2.trimmed.fastq.gz"
	output:
		assembly="{sample}/02_assembly/{sample}.assembled.fasta"
	params:
		assembly_pattern=".complete.graph1.1.path_sequence.fasta",
		additional_params = config["Assembly"]["GetOrganelle_AdditionalParameters"]
	log:
		"{sample}/02_assembly/{sample}_02_Assembly.log"
	resources:
		mem_mb=lambda wildcards, input: max(int(os.stat(input[-1]).st_size*5/1000000),config["Assembly"]["mem_mb"]),
		time=config["Assembly"]["time"]
	threads: 16
	conda:
		"../envs/GetOrganelle.yaml"
	shell:
		"""
		get_organelle_from_reads.py -1 {input.read1} -2 {input.read2} \\
			-o {wildcards.sample}/02_assembly -R {config[Assembly][MaxRounds]} -k {config[Assembly][KmerSpades]} -P {config[Assembly][PreGrouped]} \\
			-F {config[Assembly][OrganelleDatabasesAssembly]} --config-dir {config[Assembly][organelle_database_folder]}\\
			-t {threads} --prefix {wildcards.sample} --overwrite --reduce-reads-for-coverage 1000 {params.additional_params} 2> {log}
		scp {wildcards.sample}/02_assembly/*{params.assembly_pattern} {output.assembly}
		sed -i '1s/.*/>>{wildcards.sample}/' {output.assembly}
		"""



