import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

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


correct_annotation_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])
incorrect_annotation_output = os.path.join(snakemake.config["workdir"],snakemake.output[1])

error_log_file = os.path.join(snakemake.config["workdir"],snakemake.log[0])
sys.stderr = sys.stdout = open(error_log_file, "w+")

##### The Code

output_correct = []
output_incorrect = []
features_gene = []
translation_correct = True

for feature in annotation.features:
    if (feature.type != "source"):
        if (feature.type == "gene"):
            # Encounter new gene, add to the output unless there was an error in the CDS
            if (features_gene != [] and translation_correct):
                output_correct += features_gene
            else:
                output_incorrect += features_gene
            features_gene = []
            translation_correct = True 
            
            # Choose with annatator to use, default is blatX
            annotators = feature.qualifiers["annotator"][0].split(";")[0]
            if (len(annotators.split(",")) > 1):
                if ("blatX" in annotators):
                    annotator = "blatX"
                else:
                    annotator = annotators.split(",")[0]
            else:
                annotator = annotators
                    
        elif (feature.type in ["CDS","tRNA"]):
            transcript_seq = feature.extract(genome).seq
            translated_seq = feature.extract(genome).seq.translate()
            
            # Check if CDS regions are fine
            try:
                error_msg = "Feature type '%s' at position start '%d'." %(feature.type, feature.location.parts[0].start + 1)
                if (len(transcript_seq) % 3 != 0): 
                    raise NotDivisbleBy3Error(error_msg)
                if (len(transcript_seq)/3 != len(translated_seq)): 
                    raise TranslatedSequenceDiffersFromTranscriptSequenceError(error_msg)
                if (translated_seq[0] != "M"): 
                    raise DoesntStartWithStartCodonError(error_msg)
                if (translated_seq[-1] != "*"):
                    raise DoesntEndWithStopCodonError(error_msg)
                    
            except NotDivisbleBy3Error as ve:
                print(f"ERROR! Sequence is not divisible by 3: {ve}")
                translation_correct = False

            except TranslatedSequenceDiffersFromTranscriptSequenceError as fnfe:
                print(f"ERROR! Translated transcript sequence differs from the translated sequence: {fnfe}")
                translation_correct = False

            except DoesntStartWithStartCodonError as zde:
                print(f"ERROR! Feature does not start with a start codon: {zde}")
                translation_correct = False

            except DoesntEndWithStopCodonError as e:
                print(f"ERROR! Feature does not end with a stop codon: {e}")
                translation_correct = False

        # Choose only the features with the annotator of interest
        if (annotator in feature.qualifiers["annotator"][0]):    
            if ("info" in feature.qualifiers):
                feature.qualifiers.pop("info")
            if ("annotator" in feature.qualifiers):
                feature.qualifiers.pop("annotator")    
            features_gene.append(feature)
        
if (features_gene != [] and translation_correct):
    output_correct += features_gene
else:
    output_incorrect += features_gene


##### Writing Correct GenBank Entries
record_correct = SeqRecord(seq = genome.seq, 
               id = annotation.id, 
               name = annotation.name, 
               description = annotation.description)
record_correct.annotations["molecule_type"] = "genomic DNA"
record_correct.features.extend(output_correct)

with open(correct_annotation_output, "w+") as output_correct_file:
    SeqIO.write(record_correct, output_correct_file, "genbank")
    
    
##### Writing Incorrect GenBank Entries
record_incorrect = SeqRecord(seq = genome.seq, 
               id = annotation.id, 
               name = annotation.name, 
               description = annotation.description)
record_incorrect.annotations["molecule_type"] = "genomic DNA"
record_incorrect.features.extend(output_incorrect)

with open(incorrect_annotation_output, "w+") as output_incorrect_file:
    SeqIO.write(record_incorrect, output_incorrect_file, "genbank")

sys.stderr.close()