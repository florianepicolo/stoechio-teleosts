#! /usr/bin/env python3
# author : Picolo Floriane

import csv

## programme pour fusionner les informations biomart sur le whole genome
## la recherche biomart -> Human, V106, protein_coding

files_path = "/home/fpicolo/Desktop/Pathways/Data_animals/Data_human/"


# files_path = "/Users/florianepicolo/Downloads/"

fPDB = "PDB" # 7931 correspondances PDB
fUniprot = "UniprotKB" # 21003 correspondances UniprotKB
fGenename = "Genename" # 22050 correspondances Genename
fGenedesc = "Genedescription" # 22701 correspondances GeneDescription
fEntrezgene = "Entrezgene" # 22000 correspondances Entrezgene


dico_infos = dict() #dico avec ENS en clé, et différentes infos en value


def loop_file(fbiomart):
    dico_file = dict()
    with open(files_path+"Biomart-wg-ENS-"+fbiomart+".txt") as f:
        for r in f: 
            # gérer les différentes façons dont sont écrits les fichiers biomart
            if "," in r: 
                rlist = r.strip().split(",")
                # vérifier qu'on a bien que des "," et pas un mélange des deux dans le fichier
                if len(rlist[0]) > 15 : 
                    rlist = r.strip().split("\t")
            else: rlist = r.strip().split("\t")

            # retirer les lignes où il n'y a pas d'informations en fonction du format du fichier biomart
            if len(rlist) > 1 and rlist[1] != '':
                if rlist[0] not in dico_file.keys(): dico_file[rlist[0]] = [rlist[1]]
                else: dico_file[rlist[0]].append(rlist[1])
            # else: 
            #     if rlist[0] not in dico_file.keys(): dico_file[rlist[0]] = ['']
            #     else: dico_file[rlist[0]].append('') 
    return dico_file

def write_file(dicoinfos, listhumagenes):
    with open("g-infos-biomart.csv", "w", newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='', escapechar=' ')
        spamwriter.writerow(["Ensembl", fEntrezgene, fUniprot, fGenename, fGenedesc, fPDB])
        for humangene in listhumagenes:
            lign = str()
            for filename in [fEntrezgene, fUniprot, fGenename, fGenedesc, fPDB]:
                if humangene in dicoinfos[filename]: 
                    lign += ','.join(list(set(dicoinfos[filename][humangene])))
                else: lign += ''
                lign += "\t"
                 
            spamwriter.writerow([humangene, lign])


if __name__ == "__main__" :
    lhumangenes = list()
    for f in [fEntrezgene, fUniprot, fGenename, fGenedesc, fPDB]:
        dico_infos[f] = loop_file(f)
        lhumangenes += dico_infos[f].keys()
        # break
    
    lhumangenes = set(lhumangenes) # 22727 unique gene avec une information au moins
    write_file(dico_infos, lhumangenes)


