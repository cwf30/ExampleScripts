#Process raw in silico PCR and HMMER results for predicting identity and T3E repertoires. Use conda env qiime2-2022.11 for full compatibility with classifiers
import os
import re
import pandas as pd
import csv
from Bio import SeqIO
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2 import Artifact

# util function
def getFiles(directory):
    # return list of files in given directory
    files = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if filename[0] != ".":
            files.append(f)
    return files

def extractHits(f, evalueCutoff):
# HMMER outputs tables that need to be processed for analysis. files with no positive hits are 13 lines long. each hit is a single line. 
# hits get reported starting on line 3. using this info, extract lines detailing all hits:
    with open(f, 'r') as fp:
        lines=fp.readlines()
        NumHits= len(lines) - 10
        genes = {}
        # get and format hits to a easier to process list of [virulence factor, e-value, gene accession number, gene description]
        for i in range(3,NumHits):
            hit = lines[i].split()
            Evalue = float(hit[4])
            geneDescription = ' '.join(hit[18:])
            if  Evalue < evalueCutoff:
                VFOC = re.search("(.*)_",hit[2])[1]
                # rename tail fibers to something shorter
                if VFOC == 'Type_1_fiber_raw':
                    VFOC = 'Rbp1'
                if VFOC == 'Type_2_fiber_raw':
                    VFOC = 'Rbp2'
                if VFOC == 'Type_3_fiber_raw':
                    VFOC = 'Rbp3'

                
                # grab genome accession number from filename
                acc = re.search("(GCF\_.*?)\_",f)[1]
                formattedHit = [acc, VFOC, Evalue, hit[0], geneDescription]
                genes = bestmatch(formattedHit, genes)
        hits = list(genes.values())
    return hits

def bestmatch(queryGene, GeneList):
    #keeps a record of the current best annotation for a given protein
    #GeneList = {"WP_3400123":[acc, VFOC, Evalue, hit[0], geneDescription]}
    
    if queryGene[3] not in GeneList:
        GeneList[queryGene[3]] = queryGene
        return GeneList
    if queryGene[2] < GeneList[queryGene[3]][2]:
        GeneList[queryGene[3]] = queryGene
    return GeneList

#Step one: take FASTA files that describe all amplicons in each genome, and split them into new FASTA files for each primer set
def buildFASTAs(f, outDir):
    goodPrimers = ["rpoD_HWANG",
                   "gapA_Hwang",
                   "CTS_YAN",
                   "pgi_YAN",
                   "gyrB_Hwang",
                   "CTS_Hwang",
                   "CTS_Morris_Sarker",
                   "gyrB_YAN",
                   "rpoD_Parkinson"]
    amps = {}
    for i in range(len(f)):
        genome = re.search("(GCF\_.*?)\_", f[i])[1]
        with open(f[i]) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                primer = "_".join(record.id.split("_")[:-2])
                if primer not in goodPrimers:
                    continue
                if primer not in amps:
                    amps[primer] = f'>{genome}\n{record.seq}\n'
                else:
                    amps[primer] = amps[primer] + f'>{genome}\n{record.seq}\n'
    for key, value in amps.items():
        print(key)
        with open(f'{outDir}/{key}.fasta', 'w') as f:
            f.write(value)

#Step two: predict taxonomic ID 
def predictID(classifierDirectory, FASTAfileDirectory, outDir):
    classifiers = os.listdir(classifierDirectory)
    fasta_files = os.listdir(FASTAfileDirectory)
    classifiers = [i for i in classifiers if not i.startswith(".")]
    fasta_files = [i.split(".")[0] for i in fasta_files if not i.startswith(".")]

    final_data = {}
    for primer in fasta_files:

        classifier = f'{primer}.qza'
        # get the classifier from the classifiers folder
        classifier_artifact = Artifact.load(f'{classifierDirectory}/{classifier}')

        sequence_artifact = Artifact.import_data(
            'FeatureData[Sequence]', f'{FASTAfileDirectory}/{primer}.fasta')
        # classify the sequence
        classification = classify_sklearn(sequence_artifact, classifier_artifact, confidence=0.9)
        # get the classification from the classification object
        classification_artifact = classification.classification
        # save the classification to a dictionary
        classification_dict = classification_artifact.view(pd.DataFrame).to_dict()
        # get keys in classification_dict["Taxon"]
        amplicons = list(classification_dict["Taxon"].keys())
        parsed_data = {}
        # for each amplicon, get the taxon and confidence
        for amplicon in amplicons:
            parsed_data[amplicon] = {"confidence": classification_dict["Confidence"][amplicon],
                                        "level": []}
            # for each level in the classification, get the taxon and confidence
            taxon = classification_dict["Taxon"][amplicon]
            # create a list of numbers from 80 to 99
            ANIs = list(range(80, 100))
            # split taxon into levels
            levels = taxon.split(";")
            # split each level by "__"
            tmpDict = {}
            highest = 80
            for level in levels:
                level_split = level.split("__")
                # create a dictionary for each level with level_split[0] as the key and level_split[1] as the value
                if len(level_split) > 1:
                    tmpDict[int(level_split[0])] = level_split[1]
                    highest = int(level_split[0])

            parsed_data[amplicon]["level"].append(highest)
            parsed_data[amplicon]["level"].append(tmpDict[highest])
        final_data[primer] = parsed_data

    genomeAccessions = []
    with open(f'{outDir}/predictions.csv', 'w', newline='') as csvfile:
        fieldnames = ['source', 'primer', 'resolution', 'group']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for primer in final_data:
            for amplicon in final_data[primer]:
                if amplicon not in genomeAccessions:
                    genomeAccessions.append(amplicon)
                writer.writerow({'source': amplicon, 
                'primer': primer, 
                'resolution': final_data[primer][amplicon]["level"][0],
                'group': final_data[primer][amplicon]["level"][1]})
    return(pd.DataFrame({'source': genomeAccessions}))

#Step three: build T3E repertoires 
def getEffectorRepertoires(HMMERfiles, Eval, genomeAcc, outdir):
    df = genomeAcc
    for i in range(len(HMMERfiles)): 
        hits = extractHits(HMMERfiles[i], Eval)
         # check if already a column for this virulence factor. if not, add it
        for hit in hits:
            if not hit[1] in df.columns:
                df[hit[1]] = 0
            # increment count for VFOC for each one found in that file
            df.loc[df.source == hit[0], hit[1]] += 1
    outfile = f'{outdir}/repertoires.csv'
    df.to_csv(outfile)

#Step One
FASTA_in_dir = 'PCR_FASTA'
fastaFiles = getFiles(FASTA_in_dir)
FASTA_out_dir = 'PCR_Primer_FASTA'
buildFASTAs(fastaFiles, FASTA_out_dir)

#Step Two
classifierDir = 'trained_classifiers'
prediction_out_dir = 'DataForAnalysis'
genomeAccessions = predictID(classifierDir, FASTA_out_dir, prediction_out_dir)

#Step Three
HMMER_in_dir = 'HMMER_tables'
HMMERtables = getFiles(HMMER_in_dir)
getEffectorRepertoires(HMMERtables, 10e-20, genomeAccessions, prediction_out_dir)