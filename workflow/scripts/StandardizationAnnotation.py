import os
import copy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


##### Input/output paths 
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
input_seq_rec = SeqIO.read(annotation_file_input, 'genbank')
annotation_file_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])

lsc_gene_name = snakemake.config["lsc_gene"]
lsc_gene_dir = snakemake.config["lsc_gene_dir"]
ssc_gene_name = snakemake.config["ssc_gene"]
ssc_gene_dir = snakemake.config["ssc_gene_dir"]



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
    IRA_copy = copy.deepcopy(IRA)
    IR1 = copy.deepcopy(IRB)
    IR2 = IRA_copy
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
    print("%s not found in file...skipping" %lsc_gene_name)

new_seq += IR1start_IR1end    

### SSC 
try:
    if (ssc_gene.location.strand == ssc_gene_dir):
        new_seq += IR1end_IR2start
    else:
        new_seq += IR1end_IR2start.reverse_complement()
except NameError as error:
    print("%s not found in file...skipping" %ssc_gene_name)

### Rest of sequence
new_seq += IR2start_IR2_end
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