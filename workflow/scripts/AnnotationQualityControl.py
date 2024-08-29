import os, sys, logging
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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



##### Error Classes

class NotDivisbleBy3Error(Exception):
    pass

class TranslatedSequenceDiffersFromTranscriptSequenceError(Exception):
    pass

class DoesntStartWithStartCodonError(Exception):
    pass

class DoesntEndWithStopCodonError(Exception):
    pass

##### Input/output paths 

genome_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
genome = SeqIO.read(genome_file_input, 'fasta')

annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[1])
annotation = SeqIO.read(annotation_file_input, 'genbank')

annotation_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])

plant_genes_df = pd.read_csv(snakemake.params[0], names = ["genes"])
plant_genes = list(plant_genes_df["genes"])

#error_log_file = os.path.join(snakemake.config["workdir"],snakemake.log[0])
#sys.stderr = sys.stdout = open(error_log_file, "w+")

##### The Code

### Genome length QC
if (len(genome) < 120000):
    #print("WARNING! Length of genome is %d, which is shorter than expected. Expected range is 120kb to 160kb." %len(genome))
    loggerPP.warning("WARNING! Length of genome is %d, which is shorter than expected. Expected range is 120kb to 160kb." %len(genome))

elif (len(genome) > 160000):
    #print("WARNING! Length of genome is %d, which is longer than expected. Expected range is 120kb to 160kb." %len(genome))
    loggerPP.warning("WARNING! Length of genome is %d, which is longer than expected. Expected range is 120kb to 160kb." %len(genome))

else:
    #print("Length of genome is %d." %len(genome))
    loggerPP.debug("Length of genome is %d." %len(genome))
        
num_genes = {"pcg":0,"trna":0,"rrna":0}
output = []
annotator = ""
uniq_annot = set()
for feature in annotation.features:
    if (feature.type == "gene"):
        # Obtain gene names
        gene = feature.qualifiers["gene"][0]
        if gene in plant_genes:
            plant_genes.remove(gene)
        if ("rrn" in gene):
            num_genes["rrna"] += 1
        elif ("trn" in gene):
            num_genes["trna"] += 1
        else:
            num_genes["pcg"] += 1 #protein coding genes
        
        # Choose which annotator to use, default is Chloe
        annotators = feature.qualifiers["annotator"][0].split(";")[0]
        if ("blatX" in annotators):
            annotator = "blatx"
        else:
            annotator = annotators.split(",")[0]
       

    elif (feature.type in ["CDS"]):
        transcript_seq = feature.extract(genome).seq
        translated_seq = feature.extract(genome).seq.translate()

        # Check if CDS regions are fine
        try:
            error_msg = "Feature type '%s' from gene '%s' at position start '%d'." %(feature.type, gene, feature.location.parts[0].start + 1)
            if (len(transcript_seq) % 3 != 0): 
                raise NotDivisbleBy3Error(error_msg)
            if ("translation" in feature.qualifiers) and (feature.qualifiers["translation"][0] != translated_seq[:-1]): 
                raise TranslatedSequenceDiffersFromTranscriptSequenceError(error_msg)
            if (transcript_seq[0:3] != "ATG"):
                if ((gene == "psbC" and transcript_seq[0:3] == "GTG") or
                    (gene == "ndhD" and transcript_seq[0:3] == "ACG") or 
                    (gene == "rps19" and transcript_seq[0:3] == "GTG") or 
                    (gene == "ycf1" and transcript_seq[0:3] == "GTG")): 
                    pass
                else:
                    raise DoesntStartWithStartCodonError(f"{error_msg} Starts with '{transcript_seq[0:3]}'")
            if (translated_seq[-1] != "*"):
                raise DoesntEndWithStopCodonError(error_msg)

        except NotDivisbleBy3Error as ndb3e:
            #print(f"Sequence is not divisible by 3: {ndb3e}")
            loggerPP.warning(f"Sequence is not divisible by 3: {ndb3e}")

        except TranslatedSequenceDiffersFromTranscriptSequenceError as tsdftse:
            #print(f"Translated transcript sequence differs from the translated sequence: {tsdftse}")
            loggerPP.warning(f"Translated transcript sequence differs from the translated sequence: {tsdftse}")

        except DoesntStartWithStartCodonError as dswsce:
            #print(f"Feature does not start with a start codon: {dswsce}")
            loggerPP.warning(f"Feature does not start with a start codon: {dswsce}")

        except DoesntEndWithStopCodonError as dewsce:
            #print(f"Feature does not end with a stop codon: {dewsce}")
            loggerPP.warning(f"Feature does not end with a stop codon: {dewsce}")

    elif (feature.type in ["misc_feature","repeat_region"]):
        if ("annotator" in feature.qualifiers):
            annotators = feature.qualifiers["annotator"][0].split(";")[0]
            if ("blatX" in annotators):
                annotator = "blatx"
            else:
                annotator = annotators.split(",")[0]

    # Choose only the features with the annotator of interest
    if ("info" in feature.qualifiers):
        feature.qualifiers.pop("info")

    if ("annotator" in feature.qualifiers):
        if (annotator in feature.qualifiers["annotator"][0]): 
            feature.qualifiers.pop("annotator")    
            output.append(feature)
    else:
        output.append(feature)


for gene in plant_genes:
    #print("Gene %s was not found in annotation." %gene)
    loggerPP.info("Gene %s was not found in annotation." %gene)

total_genes = num_genes["pcg"] + num_genes["trna"] + num_genes["rrna"]
if (total_genes < 101):
    #print("Total number of genes found: %d, which is less than expected. Expected range is between 101 to 118 genes." %total_genes)
    loggerPP.info("Total number of genes found: %d, which is less than expected. Expected range is between 101 to 118 genes." %total_genes)
elif (total_genes > 118):
    #print("Total number of genes found: %d, which is greater than expected. Expected range is between 101 to 118 genes." %total_genes)
    loggerPP.info("Total number of genes found: %d, which is greater than expected. Expected range is between 101 to 118 genes." %total_genes)

else:
    #print("Total number of genes found: %d." %total_genes)
    loggerPP.debug("Total number of genes found: %d." %total_genes)

#print("'%d' protein coding genes found, '%d' tRNA genes found, '%d' rRNA genes found." %(num_genes["pcg"], num_genes["trna"], num_genes["rrna"]))
loggerPP.debug("'%d' protein coding genes found, '%d' tRNA genes found, '%d' rRNA genes found." %(num_genes["pcg"], num_genes["trna"], num_genes["rrna"]))


##### Writing GenBank Entries
record = SeqRecord(seq = genome.seq, 
               id = annotation.id, 
               name = annotation.name, 
               description = annotation.description)
record.annotations["molecule_type"] = "genomic DNA"
record.features.extend(output)

with open(annotation_output, "w+") as output_file:
    SeqIO.write(record, output_file, "genbank")
sys.stderr.close()