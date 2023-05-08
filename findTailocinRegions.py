""" The goal of this script is to get the genomic location of all tailocins we can. We also produce a file describing where each tail fiber that we know of is,
so that tailocin regions can be conneected with specific fibers without having to re-annotate the sequences.

output files: fiberLocations.json & tailocinLocations.json

"""
import json
import os
from Bio import SeqIO
import re
#trpG is very close to end of tailocin, but trpE is not always present. in B728a, there are 25,000bp between the start of trpG to the end of trpE, so we will 
#find the start of trpG, and grab 25kbp, or as many as we can in the contig, whichever is less.

"""
pseudocode:
go through each coding sequence.
if trpG in CDs:
    if location is on complement:
        start at starting point and move forward in the contig until we reach end of contig, or 25kbp
    else start at ending point and move backward in the contig until we reach end of contig, or 25kbp
    save the genome accession, contig name, coordinates, and accession & coordinated of known tail fiber in region to csv file.
    on ROAR, use CSV file to grab the full dna sequences
     from genomic fasta file and add to a large fasta file containing all putative tailocin regions
"""
trpG = "anthranilate synthase component II"
tailocin_locations = {}  #{"genome_accession":{"contigID": {"start":2, "end":640}}
fiber_locations = {} # {"genome_accession":{"contigID": [{"protein_accession": {"start":2, "end":640}}]}}

proteomes = os.listdir('/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/data/ncbi-cds-2022-11-17')
# get full path of every file in proteomes
proteomes = [os.path.join('/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/data/ncbi-cds-2022-11-17', proteome)
             for proteome in proteomes if proteome[0] != '.']

def getContig(seqID):
    contig = seqID.split('_')
    return('_'.join(contig[0:2]).split('|')[1])

def getLocation(seqDescription):
    location = seqDescription.split('location=')[1].split(']')[0]
    complement = False
    if 'complement' in location:
        location = location.split('complement(')[1][0:-1]
        complement = True
    if 'join' in location[0]:
        location = location.split('join(')[1][0:-1]
    location = location.split('..')
    start = int(re.sub('\D', '', location[0]))
    end = int(re.sub('\D', '', location[len(location)-1]))
    return(start, end, complement)

def getFiberLocations():
    with open('PROTEIN_VFOC_20.json') as f:
        protein_data = json.load(f)
    for key,value in protein_data.items():
        if "Rbp" not in value["HMMER"][0]:
            continue
        for genome in protein_data[key]["genomes"]:
            for proteome in proteomes:
                if genome not in proteome:
                    continue
                with open(proteome, 'r') as handle:
                        sequences = list(SeqIO.parse(handle, 'fasta'))
                        for i in range(len(sequences)):
                            if key not in sequences[i].id:
                                continue
                            
                            contig = getContig(sequences[i].id)
                            start, end, complement = getLocation(sequences[i].description)
                            
                            if genome not in fiber_locations:
                                fiber_locations[genome] = {contig: [{key: {"type": value["HMMER"][0],"start":start, "end":end}}]}
                            else:
                                if contig not in fiber_locations[genome]:
                                    fiber_locations[genome][contig] = [{key: {"type": value["HMMER"][0],"start":start, "end":end}}]
                                else: 
                                    fiber_locations[genome][contig].append({key: {"type": value["HMMER"][0],"start":start, "end":end}})

    with open('fiberLocations.json', 'w') as fp:
        json.dump(fiber_locations, fp, indent=4)


def getTailocinLocations():
    for proteome in proteomes:  
        with open(proteome, 'r') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
        for i in range(len(sequences)):  
            if trpG not in sequences[i].description:
                continue
            start, end, complement = getLocation(sequences[i].description)
            contig = getContig(sequences[i].id)
            nextGene = 1 if complement else -1
            EndNotFoundYet = True
            j = i+nextGene
            fragment_start = start if complement else end
            fragment_end = end if complement else start
            fragment_length = abs(start-end)
            while EndNotFoundYet:
                if j == len(sequences) or j < 0:
                    EndNotFoundYet = False
                    break
                if fragment_length > 25000:
                    EndNotFoundYet = False
                    break
                if getContig(sequences[j].id) != contig:
                    EndNotFoundYet = False
                    break

                s,e,c = getLocation(sequences[j].description)
                if complement:
                    fragment_length = e - fragment_start
                    fragment_end = e
                else:
                    fragment_length = fragment_start - s
                    fragment_end = s
                j += nextGene
                
                    
            

            genome = proteome.split('/')
            genome = genome[len(genome)-1].split('_')
            genome = '_'.join(genome[0:2])

            if genome not in tailocin_locations:
                tailocin_locations[genome] = {contig: [{"start":min([fragment_end,fragment_start]), "end":max([fragment_end,fragment_start]), "complement": complement}]}
            else:
                if contig not in tailocin_locations[genome]:
                    tailocin_locations[genome][contig] = [{"start":min([fragment_end,fragment_start]), "end":max([fragment_end,fragment_start]), "complement": complement}]
                else: 
                    tailocin_locations[genome][contig].append({"start":min([fragment_end,fragment_start]), "end":max([fragment_end,fragment_start]), "complement": complement})
    
    with open('tailocinLocations.json', 'w') as fp:
        json.dump(tailocin_locations, fp, indent=4)

getTailocinLocations()


