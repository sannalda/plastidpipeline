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

lsc_gene_name = snakemake.config["Standardization"]["lsc_gene"]
lsc_gene_dir = snakemake.config["Standardization"]["lsc_gene_dir"]
ssc_gene_name = snakemake.config["Standardization"]["ssc_gene"]
ssc_gene_dir = snakemake.config["Standardization"]["ssc_gene_dir"]



##### Searching for genes and invertible regions
for i in input_seq_rec.features:
    if (i.type == "gene"):
        if (i.qualifiers["gene"][0]) == lsc_gene_name:
            lsc_gene = i
        if (i.qualifiers["gene"][0]) == ssc_gene_name:
            ssc_gene = i
    if (i.type == "repeat_region" and "IRA" in i.qualifiers["note"][0]):
        IRA = i
    if (i.type == "repeat_region" and "IRB" in i.qualifiers["note"][0]):
        IRB = i

if (IRB.location.start < IRA.location.start):
    IR1 = copy.deepcopy(IRB)
    IR2 = copy.deepcopy(IRA)
else:
    IR1 = IRA
    IR2 = IRB
    
start_IR1start = input_seq_rec.seq[:IR1.location.start]
IR1start_IR1end = input_seq_rec.seq[IR1.location.start:IR1.location.end]
IR1end_IR2start = input_seq_rec.seq[IR1.location.end:IR2.location.start]
IR2start_IR2_end = input_seq_rec.seq[IR2.location.start:IR2.location.end]
IR2end_end = input_seq_rec.seq[IR2.location.end:]



##### Standardization (including reverse complementing)

new_seq = ""

### LSC 
try:
    if (lsc_gene.location.strand == lsc_gene_dir):
        new_seq += start_IR1start
    else:
        new_seq += start_IR1start.reverse_complement()
except NameError as error:
    loggerPP.error("%s not found in annotation file %s...skipping" %(lsc_gene_name,annotation_file_input))

### IRA shuold be first, then IRB (they are basically complements of each other)
if (IRB.location.start < IRA.location.start):
    new_seq += IR2start_IR2_end
else:
    new_seq += IR1start_IR1end

### SSC 
try:
    if (ssc_gene.location.strand == ssc_gene_dir):
        new_seq += IR1end_IR2start
    else:
        new_seq += IR1end_IR2start.reverse_complement()
except NameError as error:
    loggerPP.error("%s not found in annotation file %s...skipping" %(ssc_gene_name,annotation_file_input))

### IRB
if (IRB.location.start < IRA.location.start):
    new_seq += IR1start_IR1end
else:
    new_seq += IR2start_IR2_end

### Rest of Sequence
new_seq += IR2end_end



##### Output

record = SeqRecord(
    new_seq,
    id=snakemake.wildcards["sample"],
    name="FastaStandardized",
    description=""
)
# NOTE: The header of the FASTA file must contain no blanks/spaces (for some odd reason) i.e. ">Am09_Chloroplasts" is allowed, but not ">Am09 Chloroplasts"

with open(annotation_file_output, "w+") as result_file:
    SeqIO.write(record, result_file, "fasta")