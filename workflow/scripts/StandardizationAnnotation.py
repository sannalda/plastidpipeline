import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


##### Input/output paths 
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
input_seq_rec = SeqIO.read(annotation_file_input, 'genbank')
annotation_file_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])


##### Searching for genes and invertible regions
for i in input_seq_rec.features:
    if (i.type == "gene"):
        if (i.qualifiers["gene"][0]) == "psbA":
            psbA = i
        if ((i.qualifiers["gene"][0]) == "rrn23") and (i.location.strand == 1):
            rrn23_forward = i
        if ((i.qualifiers["gene"][0]) == "rrn23") and (i.location.strand == -1):
            rrn23_reverse = i
        if (i.qualifiers["gene"][0]) == "ccsA":
            ccsA = i
        if (i.qualifiers["gene"][0]) == "ndhF":
            ndhF = i
    if (i.type == "repeat_region" and "IRA" in i.qualifiers["note"][0]):
        IRA = i
    if (i.type == "repeat_region" and "IRB" in i.qualifiers["note"][0]):
        IRB = i

start_IRAstart = input_seq_rec.seq[:IRA.location.start]
IRAstart_IRAend = input_seq_rec.seq[IRA.location.start:IRA.location.end]
IRAend_IRBstart = input_seq_rec.seq[IRA.location.end:IRB.location.start]
IRBstart_IRBend = input_seq_rec.seq[IRB.location.start:IRB.location.end]
IRBend_end = input_seq_rec.seq[IRB.location.start:]



##### Reverse Complementing

new_seq = ""

### psbA should be -1
# assert(psbA.location.start < IRA.location.start):
try:
    if (psbA.location.strand == -1):
        new_seq += start_IRAstart
    else:
        new_seq += start_IRAstart.reverse_complement()
except NameError as error:
    print("psbA not found in file...skipping")


### rrn23 forward should be +1
# assert(IRA.location.start < rrn23_forward.location.start < IRA.location.end):
try:
    if (rrn23_forward.location.strand == 1):
        new_seq += IRAstart_IRAend
    else:
        new_seq += IRAstart_IRAend.reverse_complement()
except NameError as error:
    print("rrn23 forward not found in file...skipping")
    

### ccsA should be -1, ndhF should be +1, ccsA < ndhF
try:
    if (ccsA.location.start < ndhF.location.start):
        new_seq += IRAend_IRBstart
    else:
        new_seq += IRAend_IRBstart.reverse_complement()
except NameError as error:
    print("ccsA or ndhF not found in file...skipping")
# assert(IRA.location.end < ccsA.location.start < ndhF.location.start < IRB.location.start)


### rrn23 reverse should be +1
# assert(IRB.location.start < rrn23_reverse.location.start < IRB.location.end)
try:
    if (rrn23_reverse.location.strand == -1):
        new_seq += IRBstart_IRBend
    else:
        new_seq += IRBstart_IRBend.reverse_complement()
except NameError as error:
    print("rr23 reverse not found in file...skipping")

new_seq += IRBend_end

record = SeqRecord(
    new_seq,
    id="Am09",
    name="FastaStandardized",
    description=""
)
# NOTE: The header of the FASTA file must contain no blanks/spaces (for some odd reason) i.e. ">Am09_Chloroplasts" is allowed, but not ">Am09 Chloroplasts"



##### Output

#annotation_filename_standardized = '/Users/SJAnnaldasula/Documents/BGBM/Plastid/%s_standardized.fasta' %os.path.splitext(annotation_filename)[0]
with open(annotation_file_output, "w+") as result_file:
    SeqIO.write(record, result_file, "fasta")