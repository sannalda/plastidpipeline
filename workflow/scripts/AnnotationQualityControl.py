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

##### Input/output paths 
genome_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[1])
genome = SeqIO.read(genome_file_input, 'genbank')

annotation_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])

plant_genes_df = pd.read_csv(snakemake.params[0], names = ["genes"])
plant_genes = list(plant_genes_df["genes"])

alt_starts = list(snakemake.params[1].split(", "))

anno_pref = snakemake.params[2]

#print(alt_starts, " ", anno_pref)

#error_log_file = os.path.join(snakemake.config["workdir"],snakemake.log[0])
#sys.stderr = sys.stdout = open(error_log_file, "w+")

##### Genome length check

loggerPP.info("Length of the genome [bp]: %d." %len(genome))
#print("Length of the genome [bp]:\t %d" %len(genome))
if (len(genome) < 120000):
    loggerPP.warning("Genome shorter than expected [120-160 kb].")
    #print("WARNING! Genome shorter than expected (range: 120-160 kb).")
if (len(genome) > 160000):
    loggerPP.warning("Genome longer than expected [120-160 kb].")
    #print("WARNING! Genome longer than expected (range: 120-160 kb).")


##### Annotation types check

loggerPP.debug("Counting annotations ...")
#print("Counting annotations ...")

unique_feature_count = []
featuretypes = set([genome.features[i].type for i in range(len(genome.features))])
for k in featuretypes:
    if k in ("gene", "CDS", "tRNA", "rRNA"):
        unique_feature_count.append(
        len(set([l.qualifiers["gene"][0] 
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    elif k in ("exon", "intron"):
        unique_feature_count.append(
        len(set([(l.qualifiers["gene"][0], l.qualifiers["number"][0])
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    elif k in ("misc_feature", "repeat_region"):
        unique_feature_count.append(
        len(set([l.qualifiers["note"][0] 
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    else:
        unique_feature_count.append(
            sum([genome.features[j].type == k for j in range(len(genome.features))]) )

unicounts = dict(sorted(dict(zip(featuretypes, unique_feature_count)).items()))

allcounts = dict(sorted(dict(zip(featuretypes, 
                [sum([genome.features[j].type == k for j in range(len(genome.features))]) 
                 for k in set([genome.features[i].type for i in range(len(genome.features))])])).items()))
#print(allcounts)

outstr = ""
for key, value in unicounts.items():
    outstr += "\t %s:\t %i \n" %(key, value)

loggerPP.info("Number of unique features found: \n %s" %outstr)
#print("Number of distinct features found: \n", outstr)

if unicounts.get("gene", 0) < 101:
    loggerPP.warning("Less genes than expected (range: 101-118).")
    #print("WARNING! Less genes than expected (range: 101-118).")
elif unicounts.get("gene", 0) > 118:
    loggerPP.warning("More genes than expected (range: 101-118).")
    #print("WARNING! More genes than expected (range: 101-118).")

if unicounts.get("misc_feature", 0) < 2:
    loggerPP.warning("Missing LSC/SSC annotation?")
    #print("WARNING! Missing LSC/SSC annotation?")
elif unicounts.get("misc_feature", 0) > 2:
    loggerPP.warning("Extra misc_feature (e.g. LSC, SSC) annotation?")
    #print("WARNING! Extra misc_feature (e.g. LSC, SSC) annotation?")

if unicounts.get("repeat_region", 0) < 2:
    loggerPP.warning("Missing IR annotation?")
    #print("WARNING! Missing IR annotation?")
elif unicounts.get("repeat_region", 0) > 2:
    loggerPP.warning("Extra repeat_region (e.g. IR) annotation?")
    #print("WARNING! Extra repeat_region (e.g. IR) annotation?")

if unicounts.get("source", 0) != 1:
    loggerPP.warning("Source annotation missing or corrupted?")
    #print("WARNING! Source annotation missing or corrupted?")

if unicounts.get("rRNA", 0) == 1:
    loggerPP.warning("No rRNA annotations found.")
    #print("WARNING! No rRNA annotations found.")
    
if unicounts.get("tRNA", 0) == 1:
    loggerPP.warning("No tRNA annotations found.")
    #print("WARNING! No tRNA annotations found.")

##### Genes found or missing

genes_found = set([k.qualifiers["gene"][0] 
                   for k in [genome.features[i] for i in range(len(genome.features)) 
                             if genome.features[i].type =="gene"]])
outstr_f = ""
outstr_x = ""
outstr_m = ""
for gene in genes_found:
    if gene in plant_genes:
        plant_genes.remove(gene)
        outstr_f += gene + ", "
    else:
        outstr_x += gene + ", "

for gene in plant_genes:
    outstr_m += gene + ", "

loggerPP.info("Genes found: %s" %outstr_f[:-2])
#print("Genes found: ", outstr_f[:-2])

loggerPP.info("Genes not found: %s" %outstr_m[:-2])
#print("Genes not found: ", outstr_m[:-2])

if len(outstr_x) >0:
    loggerPP.warning("Please verify these gene names: %s" %outstr_x[:-2])
    #print("WARNING! Please verify these gene names: ", outstr_x[:-2])
    
##### Resolving alternative annotations & swapped IR names issue


if anno_pref in ("blat", "Chloe"):
    loggerPP.debug("Resolving alternative annotations ...")
    #print("Resolving alternative annotations ...")
    
    dupcounts = dict(sorted(dict(zip(featuretypes, [0]*len(featuretypes))).items()))
    keep_feature = [1]*len(genome.features)
    
    for k in featuretypes:    
        
        #### find indexes off all annotations of type k
        features_indexes = []
        
        ## Gene and CDS:
        if k in ("gene", "CDS", "tRNA", "rRNA"):
            features_found = set(i.qualifiers["gene"][0] for i in genome.features if i.type==k)
            for i in features_found:
                features_indexes.append([x for x,v in 
                                         enumerate(genome.features) 
                                         if v.type==k and v.qualifiers["gene"][0]==i])
        ## Exon & Intron:   
        elif k in ("exon", "intron"):
            features_found = set((i.qualifiers["gene"][0], i.qualifiers["number"][0]) 
                                 for i in genome.features if i.type==k)
            for i in list(features_found):
                features_indexes.append([x for x,v in 
                                         enumerate(genome.features) 
                                         if v.type==k
                                         and v.qualifiers["gene"][0]==i[0]
                                         and v.qualifiers["number"][0]==i[1]])
        ## LSC/SSC and IR:
        elif k in ("misc_feature", "repeat_region"):
            features_found = set(i.qualifiers["note"][0] for i in genome.features if i.type==k)
            for i in features_found:
                features_indexes.append([x for x,v in 
                                         enumerate(genome.features) 
                                         if v.type==k and v.qualifiers["note"][0]==i])
        ## Source or other
        else:
            features_found = k
            features_indexes.append([x for x,v in enumerate(genome.features) if v.type==k])
            
            
        featix = dict(zip(features_found, features_indexes))
        
        #### find and compare alternatives
        
        
        IRAnote = ""
        IRAswitch = False
        
        for feat, ixs in featix.items():
    
            
            ## if there is more than one annotation with the same feature name:
            if len(ixs) > 1:
                for i,j in zip(ixs[:-1], ixs[1:]):
                    # for each pair where list distance is smaller than 6 & start less than 1 kb apart:
                    if (j-i < 6 and (genome.features[j].location.start 
                                    - genome.features[j].location.start) < 1000):
                        pair = True
                            
                        ## QC comparison for CDS annotations
                        if (k == "CDS" and genome.features[i].qualifiers.get("translation",[0])[0] != 0):
                            
                            comparevec = [0]*8
                            
                            for fi in [0,1]:
                                
                                feat_seq = genome.features[(i,j)[fi]].extract(genome).seq
                                feat_trans = genome.features[(i,j)[fi]].extract(genome).seq.translate()
                                
                                # starts with a start codon
                                if feat_seq[0:3] != "ATG":
                                    if feat_seq[0:3] in alt_starts:
                                        comparevec[4*fi] +=.5
                                    else:
                                        comparevec[4*fi] +=1
                                    
                                # length divisible by 3
                                if (len(feat_seq) % 3 != 0):
                                    comparevec[4*fi+1] +=1
                                    
                                # transcript corresponds with sequence
                                if (genome.features[i].qualifiers.get("translation",[0])[0] 
                                    != feat_trans[:-1]):
                                    comparevec[4*fi+2] +=1
                                    
                                # ends with a stop codon
                                if feat_trans[-1] != "*":
                                    comparevec[4*fi+3] +=1
                                
                            # compare if one of the alternatives contains a problem & keep the other
                            if sum(comparevec[0:3]) > sum(comparevec[4:8]):
                                keep_feature[i] = 0
                                dupcounts[k] += 1
                                
                                loggerPP.warning("Removed alternative annotation for CDS %s (quality)." %feat)
                                #print("WARNING! Removed alternative annotation for CDS %s (quality)." %(feat))
                                 
                            elif sum(comparevec[0:3]) < sum(comparevec[4:8]):
                                keep_feature[j] = 0
                                dupcounts[k] += 1
                                
                                loggerPP.warning("Removed alternative annotation for CDS %s (quality)." %feat)
                                #print("WARNING! Removed alternative annotation for CDS %s (quality)." %(feat))
                                
                            else:
                                # preferred annotator in case of a tie
                                if anno_pref in genome.features[i].qualifiers["annotator"]:
                                   keep_feature[j] = 0 
                                   dupcounts[k] += 1
                                else:
                                   keep_feature[i] = 0 
                                   dupcounts[k] += 1
                                    
                                loggerPP.warning("Removed alternative annotation for CDS %s (annotator)." %feat)
                                #print("WARNING! Removed alternative annotation for CDS %s (annotator)." %(feat))

                                    
                        # for non-CDS / untranslatable annotations
                        else:
                            # preferred annotator
                            if anno_pref in genome.features[i].qualifiers["annotator"]:
                               keep_feature[j] == 0 
                               dupcounts[k] += 1
                               
                            else:
                               keep_feature[i] = 0 
                               dupcounts[k] += 1
                            
                            if k in ("exon", "intron"):
                                loggerPP.warning("Removed alternative annotation for %s %s (annotator)." 
                                                 %(k, str(feat[0])+ " " + str(feat[1])))
                                
                            else:
                                loggerPP.warning("Removed alternative annotation for %s %s (annotator)." %(k, feat))
                                #print("WARNING! Removed alternative annotation for %s %s (annotator)." %(k, feat))
                        
            #### IR naming & direction correction Note: currently only works if IR is not 2x annotated.
            
            if k == "repeat_region" and IRAswitch == False:
                IRAswitch = True
                
                IRAi = [i for i in [value for key, value in featix.items() if "IRA" in key][0] if keep_feature[i] == 1]
                IRBi = [i for i in [value for key, value in featix.items() if "IRB" in key][0] if keep_feature[i] == 1]
                
                if (len(IRAi) ==1 and len(IRBi) ==1):
                    loggerPP.debug("Standardizing IR annotations ...")
                    
                    IRAi = IRAi[0]
                    IRBi = IRBi[0]
                    if genome.features[IRAi].location.start > genome.features[IRBi].location.start:
                        noteA = genome.features[IRAi].qualifiers["note"]
                        noteB = genome.features[IRBi].qualifiers["note"]
                        genome.features[IRAi].qualifiers["note"] = noteB
                        genome.features[IRBi].qualifiers["note"] = noteA
                        
                        genome.features[IRAi].strand = -1
                        genome.features[IRBi].strand = 1
                        
                        loggerPP.info("Switched names of the inverted repeats.")
                        #print("WARNING! Switched names of the inverted repeats.")
                    else:
                        genome.features[IRAi].strand = 1
                        genome.features[IRBi].strand = -1
                        
else:

    IRAi = [i for i in range(len(genome.features)) if "IRA" in genome.features[i].qualifiers.get("note")]
    IRBi = [i for i in range(len(genome.features)) if "IRB" in genome.features[i].qualifiers.get("note")]
    
    if (len(IRAi) ==1 and len(IRBi) ==1):
        loggerPP.debug("Standardizing IR annotations ...")
        
        IRAi = IRAi[0]
        IRBi = IRBi[0]
        if genome.features[IRAi].location.start > genome.features[IRBi].location.start:
            noteA = genome.features[IRAi].qualifiers["note"]
            noteB = genome.features[IRBi].qualifiers["note"]
            genome.features[IRAi].qualifiers["note"] = noteB
            genome.features[IRBi].qualifiers["note"] = noteA
            
            genome.features[IRAi].strand = -1
            genome.features[IRBi].strand = 1
            
            loggerPP.info("Switched names of the inverted repeats.")
            #print("WARNING! Switched names of the inverted repeats.")
        else:
            genome.features[IRAi].strand = 1
            genome.features[IRBi].strand = -1


##### checking CDSes, compiling new features list & removing annotator + info qualifiers
loggerPP.debug("Checking CDS annotations ...")
#print("Checking CDS annotations ...")

output = []
for i in range(len(genome.features)):
    if keep_feature[i] == 1:
        feat = genome.features[i]
        
        if feat.type == "CDS":
            qc_out = ""
            
            feat_seq = feat.extract(genome).seq
            feat_trans = feat.extract(genome).seq.translate()
            
            # starts with a start codon
            if (feat_seq[0:3] != "ATG"):
                if feat_seq[0:3] in alt_starts:
                    qc_out += str(" - has alternative start codon %s" %feat_seq[0:3])
                else:
                    qc_out += str(" - starts with codon %s" %feat_seq[0:3])
                
            # length divisible by 3
            if (len(feat_seq) % 3 != 0):
                qc_out += " - has a length not a multiple of 3"
                
            # transcript corresponds with sequence
            if (genome.features[i].qualifiers.get("translation",[0])[0] 
                != feat_trans[:-1]):
                if genome.features[i].qualifiers.get("translation",[0])[0] == 0:
                    qc_out += " - has no transcript"
                else:
                    qc_out += " - has a transcript not corresponding to its sequence"
                
            # ends with a stop codon
            if feat_trans[-1] != "*":
                qc_out += " - does not end with a known stop codon"
            
            if len(qc_out) >0:
                loggerPP.warning("The CDS %s starting at %i%s. Please verify!" 
                                 %(feat.qualifiers["gene"][0], feat.location.start, qc_out))
                
                '''print("WARNING! The CDS %s starting at %i%s. Please verify!" 
                      %(feat.qualifiers["gene"][0], feat.location.start+1, qc_out))'''
        
        if ("annotator" in feat.qualifiers):
            feat.qualifiers.pop("annotator")
        if ("info" in feat.qualifiers):
            feat.qualifiers.pop("info")     
        output.append(feat)

##### Plausibility check    
unique_feature_count_out = []
for k in featuretypes:
    if k in ("gene", "CDS", "tRNA", "rRNA"):
        unique_feature_count_out.append(
        len(set([l.qualifiers["gene"][0] 
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    elif k in ("exon", "intron"):
        unique_feature_count_out.append(
        len(set([(l.qualifiers["gene"][0], l.qualifiers["number"][0])
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    elif k in ("misc_feature", "repeat_region"):
        unique_feature_count_out.append(
        len(set([l.qualifiers["note"][0] 
                 for l in [genome.features[j] for j in range(len(genome.features)) 
                           if genome.features[j].type ==k]])) )
    else:
        unique_feature_count_out.append(
            sum([genome.features[j].type == k for j in range(len(genome.features))]) )

unicounts2 = dict(sorted(dict(zip(featuretypes, unique_feature_count_out)).items()))

if unicounts != unicounts2:
    loggerPP.error("ERROR: Lost unique features during annotation QC – check code!")
    #print("ERROR: Lost unique features during annotation QC – check code!")
    print(unicounts)
    print(unicounts2)
        
outcounts = dict(sorted(dict(zip(featuretypes, 
                [sum([output[j].type == k for j in range(len(output))]) 
                 for k in set([output[i].type for i in range(len(output))])])).items()))
outstr = ""
for key, value in outcounts.items():
    outstr += "\t %s:\t %i \n" %(key, value)

loggerPP.info("Total number of annotations (with IR duplicates): \n %s" %outstr)
#print("Total number of annotations (with IR duplicates): \n", outstr)

for key in featuretypes:
    if not outcounts[key] + dupcounts[key] == allcounts[key]:
        loggerPP.error("ERROR: Lost %i %s(s) during annotation QC – check code!" 
                       %(allcounts[key] - outcounts[key] - dupcounts[key], key))
        print(key)
        print("in: %s " %allcounts[key])
        print("out: %s " %outcounts[key])
        print("rem: %s " %dupcounts[key])
            
        '''print("ERROR: Lost %i %s(s) during annotation QC – check code!" %(
            allcounts[key] - outcounts[key] - dupcounts[key], key))'''

##### writing new Genbank file
loggerPP.debug("Writing to file ...")
#print("Writing to file ...")

record = SeqRecord(seq = genome.seq, 
               id = genome.id, 
               name = genome.name, 
               description = genome.description)
record.annotations["molecule_type"] = "genomic DNA"
record.features.extend(output)

with open(annotation_output, "w+") as output_file:
    SeqIO.write(record, output_file, "genbank")
sys.stderr.close()

