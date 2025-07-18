import os, copy, logging

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


##### Input/output paths 
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
input_seq_rec = SeqIO.read(annotation_file_input, 'genbank')
annotation_file_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])

lsc_gene_name = snakemake.config["Standardize"]["lsc_gene"]
lsc_gene_dir = snakemake.config["Standardize"]["lsc_gene_dir"]
ssc_gene_name = snakemake.config["Standardize"]["ssc_gene"]
ssc_gene_dir = snakemake.config["Standardize"]["ssc_gene_dir"]


###### Locating the four plastome parts & the two indicator genes
std_features = ["IRA", "IRB", "LSC", "SSC", "lsc_gene", "ssc_gene"]
feature_dict = {}

for i in input_seq_rec.features:
    if (i.type == "repeat_region" and "IRA" in i.qualifiers["note"][0]):
        feature_dict["IRA"] = i
    if (i.type == "repeat_region" and "IRB" in i.qualifiers["note"][0]):
        feature_dict["IRB"] = i
    if (i.type == "misc_feature" and "LSC" in i.qualifiers["note"][0]):
        feature_dict["LSC"] = i
    if (i.type == "misc_feature" and "SSC" in i.qualifiers["note"][0]):
        feature_dict["SSC"] = i
    if (i.type == "CDS" and i.qualifiers["gene"][0]) == lsc_gene_name:
        feature_dict["lsc_gene"] = i
    if (i.type == "CDS" and i.qualifiers["gene"][0]) == ssc_gene_name:
        feature_dict["ssc_gene"] = i
    

### make sure none of them is missing
        
missing = [i for i in std_features if i not in feature_dict]

try:
    assert(len(missing) == 0)        
except AssertionError as error:
    loggerPP.error("%s not found in annotation file %s – skipping standardization" %(", ".join(missing),annotation_file_input))
    #print("%s not found in annotation file %s – skipping standardization" %(", ".join(missing),annotation_file_input))

### if IRB before IRA, switch how we address them

if feature_dict["IRB"].location.start < feature_dict["IRA"].location.start:
    IR1_seq = feature_dict["IRB"].extract(input_seq_rec).seq
    IR2_seq = feature_dict["IRA"].extract(input_seq_rec).seq
else:
    IR1_seq = feature_dict["IRA"].extract(input_seq_rec).seq
    IR2_seq = feature_dict["IRB"].extract(input_seq_rec).seq


####### Sequence standardization

new_seq = ""

### LSC & IR1
if (feature_dict["lsc_gene"].location.strand == lsc_gene_dir):
    lsc_seq = feature_dict["LSC"].extract(input_seq_rec).seq
else:    
    lsc_seq = feature_dict["LSC"].extract(input_seq_rec).seq.reverse_complement()

new_seq += lsc_seq + IR1_seq

### SSC + IR2
if (feature_dict["ssc_gene"].location.strand == ssc_gene_dir):
    ssc_seq = feature_dict["SSC"].extract(input_seq_rec).seq
else:    
    ssc_seq = feature_dict["SSC"].extract(input_seq_rec).seq.reverse_complement()

new_seq += ssc_seq + IR2_seq


####### Check total length to make sure it worked

try:
    assert(len(new_seq) == len(input_seq_rec))        
except AssertionError as error:
    loggerPP.error("Standardized sequence differs by %s nt from input length for annotation file %s" %(len(input_seq_rec)-len(new_seq),annotation_file_input))
    raise


##### Output

record = SeqRecord(
    new_seq,
    id=snakemake.wildcards["sample"],
    name="FastaStandardized",
    description=""
)
# NOTE: The header of the FASTA file must contain no blanks/spaces (for some odd reason) i.e. ">SAMPLE_Chloroplasts" is allowed, but not ">SAMPLE Chloroplasts"

with open(annotation_file_output, "w+") as result_file:
    SeqIO.write(record, result_file, "fasta")