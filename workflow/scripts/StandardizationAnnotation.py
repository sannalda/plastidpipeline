import os

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

start_IRAstart = input_seq_rec.seq[:IRA.location.start]
IRAstart_IRAend = input_seq_rec.seq[IRA.location.start:IRA.location.end]
IRAend_IRBstart = input_seq_rec.seq[IRA.location.end:IRB.location.start]
IRBstart_IRBend = input_seq_rec.seq[IRB.location.start:IRB.location.end]
IRBend_end = input_seq_rec.seq[IRB.location.end:]



##### Reverse Complementing

new_seq = ""

### LSC 
try:
    if (lsc_gene.location.strand == lsc_gene_dir):
        new_seq += start_IRAstart
    else:
        new_seq += start_IRAstart.reverse_complement()
except NameError as error:
    print("%s not found in file...skipping" %lsc_gene_name)

new_seq += IRAstart_IRAend    

### SSC 
try:
    if (ssc_gene.location.strand == ssc_gene_dir):
        new_seq += IRAend_IRBstart
    else:
        new_seq += IRAend_IRBstart.reverse_complement()
except NameError as error:
    print("%s not found in file...skipping" %ssc_gene_name)

new_seq += IRBstart_IRBend
new_seq += IRBend_end

record = SeqRecord(
    new_seq,
    id=snakemake.wildcards["sample"],
    name="FastaStandardized",
    description=""
)
# NOTE: The header of the FASTA file must contain no blanks/spaces (for some odd reason) i.e. ">Am09_Chloroplasts" is allowed, but not ">Am09 Chloroplasts"


##### Output

#annotation_filename_standardized = '/Users/SJAnnaldasula/Documents/BGBM/Plastid/%s_standardized.fasta' %os.path.splitext(annotation_filename)[0]
with open(annotation_file_output, "w+") as result_file:
    SeqIO.write(record, result_file, "fasta")