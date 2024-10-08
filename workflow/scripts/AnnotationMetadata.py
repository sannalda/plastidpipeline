from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
import Bio
import pandas as pd
from datetime import datetime
import os, sys, logging

def setup_logger(log_file_path):
    loggerPP = logging.getLogger("PlastidPipeline_%s" %snakemake.wildcards["sample"])
    loggerPP.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    console_handler.setFormatter(formatter)
    loggerPP.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    loggerPP.addHandler(file_handler)

    return loggerPP

loggerPP = setup_logger(os.path.join(snakemake.config["workdir"],snakemake.log[0]))


##### Input Files
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[1])
annotation = SeqIO.read(annotation_file_input, 'genbank')

genome_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[2])
genome = SeqIO.read(genome_file_input, 'fasta')

try:
    metadata = pd.read_csv(os.path.join(snakemake.config["workdir"],snakemake.input[0]), sep = "\t", header = None, names = ["fields","info"], index_col = 0)
except pd.errors.ParserError as inst:
    loggerPP.error("\nThere seems to be a formatting problem in the metadata file: %s. Maybe there is an extra 'tab' somewhere. Exiting the script..." %snakemake.input[0])

metadata = metadata.T
if ("SOURCE_MOD" in metadata.columns):
    source_mod = list(metadata.columns[list(metadata.columns).index("SOURCE_MOD")+1:])
    metadata = metadata.drop(columns=source_mod + ["SOURCE_MOD"])
    #metadata["SOURCE_MOD"]["info"] = source_mod
    metadata = metadata.to_dict()
else:
    source_mod = None

#####
seq = genome.seq
locus = snakemake.wildcards["sample"]
definition = "%s accession %s plastid, complete genome" %(metadata["SPECIES"]["info"],metadata["SAMPLE"]["info"])

##### 
annotations = dict()
annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
annotations["source"] = metadata["SOURCE"]["info"]
annotations["accession"] = metadata["SAMPLE"]["info"]
annotations["organism"] = metadata["SPECIES"]["info"]
annotations["taxonomy"] = [metadata["TAXONOMY"]["info"]]
annotations["comment"] = """Preliminary file - please verify!
            Annotated plastome generated by PlastidPipeline.
            Please cite the following sources in your publication:
            1- Siddarth Annaldasula, Manuel Vera Rodriguez, Adrian Casanova 
            Chiclana, Thomas Borsch, Katja Reichel (2024): 
            Title, Journal, doi.
            2- Michael Tillich,  Pascal Lehwark,  Tommaso Pellizzer,  Elena S.
            Ulbricht-Jones,  Axel Fischer, Ralph Bock, Stephan Greiner (2017):
            GeSeq - versatile and accurate annotation of organelle genomes. 
            Nucleic Acids Research. https://doi.org/10.1093/nar/gkx391
            3- Jian-Jun Jin, Wen-Bin Yu, Jun-Bo Yang, Yu Song, Claude W. 
            dePamphlis, Ting-Shuang Yi, De-Zhu Li (2020) GetOrganelle: a fast 
            and versatile toolkit for accurate de novo assembly of organelle
            genomes. Genome Biology. https://doi.org/10.1186/s13059-020-02154-5"""
annotations["molecule_type"] = "DNA"
annotations["topology"] = "circular"

reference = Reference()
reference.authors = metadata["AUTHORS"]["info"]
reference.title = "Direct Submission"
reference.journal = metadata["INSTITUTION"]["info"]
annotations["references"] = [reference]

##### Feature Source creation
features_source_dict = {"organism":metadata["SPECIES"]["info"], 
                        "organelle":"plastid:chloroplast",
                        "molecule_type":"genomic DNA",
                        "db_xref":"taxon:"+metadata["TAXREF"]["info"]}
if (source_mod != None):
    for source_mod_line in source_mod:
        source_mod_info = source_mod_line.replace("/","").split("=")
        features_source_dict[source_mod_info[0]] = source_mod_info[1].strip('\"').strip("\'")

feature_source = SeqFeature(FeatureLocation(0, len(seq)), type="source", qualifiers = features_source_dict)


##### Output File
record = SeqRecord(seq = seq, 
                   id = locus, 
                   name = locus, 
                   description = definition, 
                   annotations = annotations, 
                   features = [feature_source])
record.annotations["data_file_division"] = "PLN"

if (annotation.features[0].type == "source"):
    annotation.features.remove(annotation.features[0])

record.features.extend(annotation.features)

SeqIO.write(record, os.path.join(snakemake.config["workdir"],snakemake.output[0]), "genbank")