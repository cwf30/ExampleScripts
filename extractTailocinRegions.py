""" This script completes the extraction of putative tailocin regions (25kbp starting at trpG), a companion script to findTailocinRegions.py"""
import json
import os
from Bio import SeqIO

with open('tailocinLocations.json') as f:
    tailocinLocations = json.load(f)
with open('fiberLocations.json') as f:
    fiberLocations = json.load(f)

def getFibersinTailocinSeq(FiberContigs, accession, contig, start, end):
    #return list of fibers that are found in a given tailocin sequence
    fibers = []
    for fiberContig in FiberContigs:
        if fiberContig in contig:
            for fiberObj in fiberLocations[accession][fiberContig]:
                for fiberName in fiberObj.keys():
                    fiber = fiberObj[fiberName]
                    if (int(fiber["start"]) < end and int(fiber["start"]) > start) or (int(fiber["end"]) < end and int(fiber["end"]) > start):
                        fibers.append(fiberName)
    return(fibers)

def reverseComplement(seq):
    return(seq.replace("A", "T").replace("C", "G").replace("T", "A").replace("G", "C")[::-1])

def extractPutativeTailocins():
    genomes = os.listdir('ncbi-genomes-2022-12-16')
    # get full path of every file in proteomes
    genomes = [os.path.join('ncbi-genomes-2022-12-16', genomes)
                for genomes in genomes if genomes[0] != '.']
    Tailocins = []
    for accession in tailocinLocations.keys():
        for genome in genomes:
            if accession not in genome:
                continue

            with open(genome, 'r') as handle:
                sequences = list(SeqIO.parse(handle, 'fasta'))
            for i in range(len(sequences)):
                for contig in tailocinLocations[accession].keys():
                    if contig not in sequences[i].id:
                        continue
                    for possibleTailocin in tailocinLocations[accession][contig]:
                        start = int(possibleTailocin["start"])
                        end = int(possibleTailocin["end"])

                        FiberContigs = []
                        if accession in fiberLocations.keys():
                            FiberContigs = fiberLocations[accession].keys()
                        fibersInTailocinSeq = getFibersinTailocinSeq(FiberContigs, accession, contig, start, end)
                        
                        tailocinSeq = sequences[i].seq[start-1:end]
                        if possibleTailocin["complement"] == "true":
                            pass
                            #tailocinSeq = reverseComplement(tailocinSeq)
                        Tailocins.append({"genome": accession,"contig":contig,"coords":f'{start}_{end}',"sequence":tailocinSeq, "fibers":("_").join(fibersInTailocinSeq)})
        
    with open("TailocinRegions.fasta", 'w') as outfile:
        for tailocin in Tailocins:
            outfile.write(f'>{tailocin["genome"]}|{tailocin["contig"]}_{tailocin["coords"]}|{tailocin["fibers"]}\n {tailocin["sequence"]}\n')

extractPutativeTailocins()