from Bio import SeqIO
import os
import json
import csv

proteomes = os.listdir('/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/data/ncbi-cds-2022-11-17')
# get full path of every file in proteomes
proteomes = [os.path.join('/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/data/ncbi-cds-2022-11-17', proteome)
             for proteome in proteomes if proteome[0] != '.']

hmmerHits = '/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/PROTEIN_VFOC_20.json'

with open(hmmerHits, 'r') as f:
    hits = json.load(f)
    fibers = {}
    tailocins = {}
    for key, value in hits.items():
        if "Rbp" not in value["HMMER"][0]:
            continue
        fibers[key] = {"type": value["HMMER"][0], "genomes": value["genomes"], "prophages":[], "seq":"", "copyNumber": 0}
        tailocins[key] = {"type": value["HMMER"][0], "genomic_regions": []}
        #print(f'{len(fibers[key]["genomes"])} genomes for fiber {key}')
        count = 0
        for genome in fibers[key]["genomes"]:
            for proteome in proteomes:
                if genome not in proteome:
                    continue
                with open(proteome, 'r') as handle:
                    sequences = list(SeqIO.parse(handle, 'fasta'))
                    for i in range(len(sequences)):
                        if key not in sequences[i].id:
                            continue
                        count += 1
                        fibers[key]["seq"] = sequences[i].seq
                        new_region = {"genome":genome}
                        gene_number = 1
                        tailocins[key]["genomic_regions"].append(new_region)
                        for j in range(i-25, i+25):
                            if len(sequences)-1 > j:
                                new_region[gene_number] = {"gene": str(sequences[j].description),"seq":str(sequences[j].seq)}
                                gene_number += 1
                        for j in range(i-25, i+25):
                            if len(sequences)-1 > j:
                                if 'capsid' in sequences[j].description:
                                    fibers[key]["prophages"].append(genome)
                                    #print('prophage ' + key)
                                    break
                        
        fibers[key]["copyNumber"] = count
        
    with open("/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/uniqueFibers.fasta", 'w') as outfile:
        for key in fibers.keys():
            if len(fibers[key]["seq"]) > 5:
                outfile.write(f'>{key} {fibers[key]["type"]}\n {fibers[key]["seq"]} \n')
    
    with open("/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/uniqueFibersType1.fasta", 'w') as outfile:
        for key in fibers.keys():
            if len(fibers[key]["seq"]) > 5 and fibers[key]["type"] == "Rbp1":
                outfile.write(f'>{key} {fibers[key]["type"]}\n {fibers[key]["seq"]} \n')
    
    with open("/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/uniqueFibersType2.fasta", 'w') as outfile:
        for key in fibers.keys():
            if len(fibers[key]["seq"]) > 5 and fibers[key]["type"] == "Rbp2":
                outfile.write(f'>{key} {fibers[key]["type"]}\n {fibers[key]["seq"]} \n')

    with open("/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/uniqueFibersType3.fasta", 'w') as outfile:
        for key in fibers.keys():
            if len(fibers[key]["seq"]) > 5 and fibers[key]["type"] == "Rbp3":
                outfile.write(f'>{key} {fibers[key]["type"]}\n {fibers[key]["seq"]} \n')

    with open("/Users/cwf30/Desktop/Code/Dissertation/Ch2_tail_fiber_analysis/FiberRegions.fasta", 'w') as outfile:
        for key in tailocins.keys():
            for fiber in tailocins[key]["genomic_regions"]:
                seq = ""
                for fiberkey in fiber:
                    if fiberkey != "genome":
                        seq = seq + fiber[fiberkey]["seq"]
                outfile.write(f'>{key}_from_{fiber["genome"]} [fiber_type={fibers[key]["type"]}] [seq_length={len(seq)}]\n {seq} \n')

    with open('fiber_regions.json', 'w') as fp:
        json.dump(tailocins, fp, indent=4)


    with open('fibersMetadata.csv', 'w', newline='') as csvfile:
        fieldnames = ['name', 'type', 'count', 'copies', 'prophage']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for key in fibers.keys():
            writer.writerow({'name': key, 'type': fibers[key]["type"], 
            'count':len(fibers[key]["genomes"]), 'copies': fibers[key]["copyNumber"],'prophage':len(fibers[key]["prophages"])/len(fibers[key]["genomes"])})
